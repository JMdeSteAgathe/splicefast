"""Interface en ligne de commande de SpliceFast."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np

from . import __version__
from .bam import DEFAULT_EXCLUDE, bam_refs, ref_offsets
from .count import (count_all, discover_samples, read_sample_sheet,
                    write_sample_sheet)
from .gtf import Annotation, parse_gtf
from .util import (Workdir, check_space, die, free_space, human, log,
                   read_json, safe_name, set_quiet, warn, winpath, write_json)


# --------------------------------------------------------------------------- #
def _samples_from_args(a) -> List[Tuple[str, str]]:
    samples: List[Tuple[str, str]] = []
    if a.samples:
        samples += read_sample_sheet(a.samples)
    if a.run:
        samples += discover_samples(a.run, recursive=not a.no_recursive,
                                    pattern=a.pattern)
    if a.bam:
        from .count import clean_sample_id
        for b in a.bam:
            samples.append((clean_sample_id(Path(b)), b))
    # dedupliquer en gardant l'ordre
    seen = set()
    out = []
    for sid, p in samples:
        if sid in seen:
            continue
        seen.add(sid)
        out.append((sid, p))
    return out


def _project(wd: Workdir, samples: Sequence[Tuple[str, str]], gtf: Optional[str],
             force_ann: bool = False, biotypes: Optional[Sequence[str]] = None):
    """Prepare references + annotation + liste globale de sites."""
    pj = wd.store / "project.json"
    if not samples:
        die("aucun echantillon : utiliser --run / --samples / --bam")
    ref_names, ref_lengths = bam_refs(samples[0][1])
    log(f"References du BAM : {len(ref_names)} sequences "
        f"(ex: {', '.join(ref_names[:4])}...)")
    ref_off = ref_offsets(ref_lengths)

    ann_path = wd.store / "annotation.npz"
    sites_path = wd.store / "sites.npy"
    if gtf and (force_ann or not ann_path.exists() or not sites_path.exists()):
        ann = parse_gtf(gtf, ref_names, ref_off, keep_biotypes=biotypes)
        ann.save(ann_path)
        np.save(winpath(sites_path), ann.site_pos.astype(np.int64))
    elif ann_path.exists():
        ann = Annotation.load(ann_path)
        if not sites_path.exists():
            np.save(winpath(sites_path), ann.site_pos.astype(np.int64))
    else:
        warn("pas de GTF fourni : pas de detection de retention d'intron "
             "(theta) ni de classification des evenements")
        ann = Annotation.empty()
        np.save(winpath(sites_path), np.zeros(0, dtype=np.int64))

    sp = np.load(winpath(sites_path))
    h = hashlib.md5(np.ascontiguousarray(sp).tobytes()).hexdigest()[:16]
    write_json(pj, dict(version=__version__, ref_names=ref_names,
                        ref_lengths=[int(x) for x in ref_lengths],
                        gtf=gtf, n_sites=int(len(sp)), sites_hash=h,
                        samples=[s for s, _ in samples],
                        bams=[b for _, b in samples]))
    return ann, ref_names, ref_lengths, ref_off, str(sites_path), h


def _load_project(wd: Workdir):
    pj = wd.store / "project.json"
    if not pj.exists():
        die(f"projet introuvable dans {wd.store} : lancer d'abord 'count'")
    meta = json.load(open(winpath(pj), encoding="utf-8"))
    ann_path = wd.store / "annotation.npz"
    ann = Annotation.load(ann_path) if ann_path.exists() else Annotation.empty()
    return meta, ann


# --------------------------------------------------------------------------- #
def cmd_discover(a) -> None:
    samples = _samples_from_args(a)
    if not samples:
        die("aucun BAM trouve")
    out = Path(a.out) if a.out else None
    log(f"{len(samples)} BAM trouves")
    tot = 0
    for sid, p in samples:
        try:
            sz = os.path.getsize(winpath(p))
        except OSError:
            sz = 0
        tot += sz
        bai = Path(str(p) + ".bai").exists() or Path(p).with_suffix(".bai").exists()
        print(f"{sid}\t{p}\t{human(sz)}\t{'bai:oui' if bai else 'bai:NON'}")
    print(f"# total : {len(samples)} echantillons, {human(tot)}", file=sys.stderr)
    if out:
        write_sample_sheet(out, samples)
        log(f"Feuille d'echantillons ecrite : {out}")


def _apply_io_opts(a) -> None:
    """Propage les reglages d'E/S, y compris aux processus fils (spawn
    Windows recharge les modules : on passe donc par l'environnement)."""
    from . import util
    if getattr(a, "io_retries", None):
        util.IO_RETRIES = int(a.io_retries)
        os.environ["SPLICEFAST_IO_RETRIES"] = str(int(a.io_retries))
    if getattr(a, "io_wait", None):
        util.IO_WAIT = float(a.io_wait)
        os.environ["SPLICEFAST_IO_WAIT"] = str(float(a.io_wait))


def cmd_count(a) -> None:
    _apply_io_opts(a)
    wd = Workdir(a.workdir).create()
    samples = _samples_from_args(a)
    ann, ref_names, ref_lengths, ref_off, sites_path, h = _project(
        wd, samples, a.gtf, force_ann=a.force_annotation,
        biotypes=a.biotypes.split(",") if a.biotypes else None)
    write_sample_sheet(wd.root / "samples.tsv", samples)
    check_space(wd.root, 20 << 30)
    opts = dict(anchor=a.anchor, min_mapq=a.min_mapq,
                exclude_flags=DEFAULT_EXCLUDE if not a.keep_duplicates
                else (DEFAULT_EXCLUDE & ~0x400),
                min_intron=a.min_intron, max_intron=a.max_intron,
                min_anchor_split=a.min_anchor_split)
    tg = None
    if getattr(a, "retention", False) and len(ann.intron_start):
        tg = _cov_targets(a, ann, ref_names, ref_lengths)
    count_all(wd, samples, ref_off, sites_path if len(np.load(winpath(sites_path))) else None,
              h, opts, threads=a.threads, force=a.force,
              min_chunk_mb=a.chunk_mb, engine=a.engine, jobs=a.jobs,
              ref_names=ref_names, cov_targets=tg)
    missing = [sid for sid, _ in samples if not wd.sample_counts(sid).exists()]
    if missing:
        warn(f"Comptage incomplet : {len(missing)} echantillon(s) sans "
             f"comptage ({', '.join(missing[:8])}"
             f"{' ...' if len(missing) > 8 else ''})")
    else:
        log("Comptage termine.")
    return missing


def _genes_from_args(a) -> Optional[List[str]]:
    genes = list(getattr(a, "genes", None) or [])
    if getattr(a, "genes_file", None):
        for line in open(winpath(a.genes_file), encoding="utf-8"):
            genes += [x for x in line.replace(",", " ").split()]
    return sorted(set(genes)) if genes else None


def _cov_targets(a, ann, ref_names, ref_lengths, genes=None):
    from .coverage import build_cov_targets
    return build_cov_targets(ann, ref_names, ref_lengths,
                             genes=genes if genes is not None else _genes_from_args(a),
                             edge=a.edge, min_len=a.min_intron_len,
                             border_width=a.border_width)


def cmd_build(a) -> None:
    from .store import build_store

    wd = Workdir(a.workdir).create()
    meta, ann = _load_project(wd)
    samples = read_sample_sheet(str(wd.root / "samples.tsv"))
    if a.samples:
        samples = read_sample_sheet(a.samples)
    missing = [sid for sid, _ in samples if not wd.sample_counts(sid).exists()]
    if missing:
        die(f"{len(missing)} echantillon(s) de la feuille sans comptage "
            f"({', '.join(missing[:5])}{' ...' if len(missing) > 5 else ''}) : "
            f"les matrices auraient des lignes vides, ce qui fausse la "
            f"reference de toute la cohorte. Relancez 'count', ou retirez-les "
            f"de la feuille d'echantillons.", code=1)
    build_store(wd, samples, ann, meta["ref_names"], meta["ref_lengths"],
                min_total=a.min_total, min_samples=a.min_samples,
                block_size=a.block_size)
    _apply_reference(a, wd)
    log("Store construit.")


def _reference_lists(a):
    """(--test-only + fichier, --reference-file) -> listes de noms."""
    names = list(getattr(a, "test_only", None) or [])
    f = getattr(a, "test_only_file", None)
    if f:
        names += [x.split("\t")[0].strip() for x in
                  open(winpath(f), encoding="utf-8").read().splitlines()
                  if x.strip() and not x.startswith("#")]
    refonly = None
    f = getattr(a, "reference_file", None)
    if f:
        refonly = [x.split("\t")[0].strip() for x in
                   open(winpath(f), encoding="utf-8").read().splitlines()
                   if x.strip() and not x.startswith("#")]
    return names, refonly


def _apply_reference(a, wd) -> bool:
    """Definit les echantillons de reference si l'utilisateur l'a demande."""
    from .store import set_reference

    names, refonly = _reference_lists(a)
    mr = int(getattr(a, "min_reads", 0) or 0)
    match = getattr(a, "reference_match", None)
    if not (names or refonly or mr or match):
        return False
    set_reference(wd, min_reads=mr, test_only=names, reference_only=refonly,
                  match=match)
    return True


def cmd_probe(a) -> None:
    """Relit une region dans le BAM et la compare au comptage stocke."""
    from .probe import probe

    wd = Workdir(a.workdir)
    try:
        _, ann = _load_project(wd)
    except Exception:  # noqa: BLE001
        ann = None
    probe(wd, ann, a.sample, a.region, bam=a.bam, min_mapq=a.min_mapq,
          anchor=a.anchor)


def cmd_why(a) -> None:
    """Pourquoi tel evenement est-il (ou n'est-il pas) dans le rapport ?"""
    from .diagnose import explain, print_explanation

    wd = Workdir(a.workdir)
    _, ann = _load_project(wd)
    rows = explain(wd, ann, a.gene, a.sample, chrom=a.chrom, start=a.start,
                   end=a.end, padj_cutoff=a.padj, delta_cutoff=a.delta,
                   nominal_p=a.nominal_p, top_per_sample=a.top_per_sample,
                   top_n=a.top)
    print_explanation(rows)
    if a.gene:
        from .coverage import gene_coverage_table
        cov = [r for r in gene_coverage_table(wd, ann, a.gene)
               if r["sample"] == a.sample]
        if cov:
            log("")
            log("module couverture (retention d'intron) pour cet echantillon :")
            for r in cov:
                log(f"  {r['chrom']}:{r['start']}-{r['end']} "
                    f"intron {r['prof_intron']:>8} / bordures {r['prof_bordures']:>8} "
                    f"ratio {str(r['ratio']):>8} (median cohorte "
                    f"{str(r['ratio_median_cohorte']):>8}, z {str(r['z_robuste']):>6}, "
                    f"padj_gene {str(r['padj_gene'])[:9]:>9}) "
                    f"exploitable={r['exploitable']}")


def cmd_reference(a) -> None:
    """Choisir qui definit la norme, sans rien recompter."""
    from .store import read_reference, set_reference

    wd = Workdir(a.workdir)
    meta = read_json(wd.store / "store.json")
    if not meta:
        die("store introuvable : lancer d'abord 'build'")
    samples = list(meta["samples"])
    reads = {q["sample"]: int(q["reads"]) for q in meta.get("qc", [])}
    if _apply_reference(a, wd):
        log("Relancer 'analyze' (et 'report') pour recalculer avec cette "
            "reference ; les comptages et la couverture sont deja faits.")
    ref = read_reference(wd, samples)
    log(f"{int(ref.sum())}/{len(samples)} echantillon(s) de reference")
    if a.list:
        for sid, r in zip(samples, ref):
            log(f"  {'reference   ' if r else 'hors ref.   '}{sid:45s} "
                f"{reads.get(sid, 0):>14,} reads")


def cmd_analyze(a) -> None:
    from .analyze import analyze

    wd = Workdir(a.workdir).create()
    cov = None
    if a.covariates:
        cov = _read_covariates(a.covariates, wd)
    analyze(wd, metrics=tuple(a.metrics.split(",")),
            n_factors=(None if a.n_factors in (None, "auto") else int(a.n_factors)),
            sub_rows=a.sub_rows, min_k_one=a.min_k_one,
            min_n_quantile=a.min_n_quantile, min_cov=a.min_cov,
            covariates=cov, robust_iter=a.robust_iter, screen=a.screen)
    log("Analyse terminee.")


def _read_covariates(path: str, wd: Workdir) -> Optional[np.ndarray]:
    import csv

    from .store import Store

    st = Store(wd)
    rows = {}
    with open(winpath(path), newline="", encoding="utf-8") as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = next(rd)
        for r in rd:
            if not r:
                continue
            rows[r[0]] = r[1:]
    cols = head[1:]
    mat = []
    for c_i, cname in enumerate(cols):
        vals = [rows.get(s, [""] * len(cols))[c_i] for s in st.samples]
        try:
            v = np.array([float(x) for x in vals])
            mat.append((v - v.mean()) / (v.std() + 1e-9))
        except ValueError:
            for lev in sorted(set(vals))[1:]:
                mat.append(np.array([1.0 if x == lev else 0.0 for x in vals]))
    if not mat:
        return None
    log(f"{len(mat)} covariable(s) ajoutee(s) au modele")
    return np.stack(mat, axis=1)


def cmd_report(a) -> None:
    from .report import make_reports

    wd = Workdir(a.workdir).create()
    _, ann = _load_project(wd)
    make_reports(wd, ann, padj_cutoff=a.padj, delta_cutoff=a.delta,
                 top_n=a.top, samples=a.sample or None,
                 shift_window=a.shift_window,
                 pseudo_exon_max=a.pseudo_exon_max,
                 panel_files=a.panel or None,
                 max_html_rows=a.html_rows,
                 cov_min_border=a.cov_min_border,
                 cov_min_ratio=a.cov_min_ratio,
                 cov_min_samples=a.cov_min_samples,
                 top_per_sample=a.top_per_sample, nominal_p=a.nominal_p,
                 cov_top_per_sample=a.cov_top,
                 harmonise=not a.no_harmonise,
                 panelapp_cache=a.panelapp, recurrent_min=a.recurrent,
                 write_excel=not a.no_excel)


def cmd_panelapp(a) -> None:
    """Telecharge (ou importe) PanelApp Australie : tous les panels."""
    from .panelapp import default_cache, download, from_json

    wd = Workdir(a.workdir).create()
    out = Path(a.out) if a.out else default_cache(wd)
    if a.from_json:
        from_json(Path(a.from_json), out)
    else:
        download(out, api=a.api, insecure=a.insecure,
                 save_json=Path(a.save_json) if a.save_json else None)
    log("Relancer `report` pour annoter les resultats avec ce cache.")


def cmd_gene(a) -> None:
    """Tableau couverture d'un gene pour tous les echantillons (sans FDR)."""
    from .coverage import gene_coverage_table

    wd = Workdir(a.workdir).create()
    _, ann = _load_project(wd)
    rows = gene_coverage_table(wd, ann, a.gene, min_border=a.min_border)
    if not rows:
        return
    out = wd.results / f"gene_{safe_name(a.gene)}.tsv"
    cols = list(rows[0].keys())
    wd.results.mkdir(parents=True, exist_ok=True)
    with open(winpath(out), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    rows.sort(key=lambda r: (r["start"], r["exploitable"] != "oui",
                             -(r["ratio"] if r["ratio"] != "" else -1)))
    cur = None
    for r in rows:
        if r["start"] != cur:
            cur = r["start"]
            log(f"--- intron {r['chrom']}:{r['start']}-{r['end']} "
                f"(mediane cohorte {r['ratio_median_cohorte']}) ---")
            shown = 0
        if (shown < a.top and r["exploitable"] == "oui") \
                or (a.sample and r["sample"] == a.sample):
            mark = "  <<<" if a.sample and r["sample"] == a.sample else ""
            if r["exploitable"] == "oui":
                log(f"  {r['sample']:34s} intron={r['prof_intron']:8.2f} "
                    f"bordures={r['prof_bordures']:8.2f} "
                    f"ratio={float(r['ratio']):.3f} "
                    f"z={float(r['z_robuste']):7.2f} "
                    f"rang={r['rang']:3d}/{r['n_exploitables']} "
                    f"p={float(r['p_value']):.2e} "
                    f"padj_gene={float(r['padj_gene']):.2e}{mark}")
            else:
                log(f"  {r['sample']:34s} intron={r['prof_intron']:8.2f} "
                    f"bordures={r['prof_bordures']:8.2f}   "
                    f"(couverture insuffisante, non classe){mark}")
            shown += 1
    log(f"tableau complet : {out}")


def cmd_retention(a) -> None:
    """Mesure la couverture des introns (module dedie, lecture ciblee)."""
    from .bam import bam_refs
    from .coverage import measure_coverage, save_coverage
    from .panels import load_panels

    wd = Workdir(a.workdir).create()
    _, ann = _load_project(wd)
    # On s'aligne sur les echantillons DU STORE : ce sont ceux qui ont ete
    # comptes et analyses. samples.tsv peut en contenir davantage (BAM
    # decouverts mais non comptes), ce qui produirait un rapport incoherent.
    samples = []
    try:
        meta = read_json(wd.store / "store.json")
        samples = list(zip(meta.get("samples", []), meta.get("bams", [])))
    except Exception:
        pass
    sheet = wd.root / "samples.tsv"
    if not samples and sheet.exists():
        samples = read_sample_sheet(str(sheet))
    elif samples and sheet.exists():
        n_sheet = len(read_sample_sheet(str(sheet)))
        if n_sheet != len(samples):
            warn(f"samples.tsv liste {n_sheet} echantillons mais le store en "
                 f"contient {len(samples)} : on utilise ceux du store")
    if not samples:
        die("aucun echantillon : lancer 'discover' puis 'count' d'abord")
    ref_names, ref_lengths = bam_refs(samples[0][1])

    genes = _genes_from_args(a) or []
    if a.panel:
        ps = load_panels(a.panel)
        if ps:
            genes += list(ps.genes.keys())
    genes = sorted(set(genes)) if genes else None
    tg = _cov_targets(a, ann, ref_names, ref_lengths, genes=genes)
    S, M = measure_coverage(wd, samples, tg, ref_lengths, threads=a.threads,
                            chunk_mb=a.chunk_mb, jobs=a.jobs, engine=a.engine,
                            ref_names=ref_names)
    save_coverage(wd, tg, samples, S, M)
    log("termine : relancer 'report' pour integrer ces resultats")


def _store_ready(wd: Workdir, samples) -> bool:
    """Le store est-il deja construit, complet et coherent avec les echantillons ?"""
    try:
        meta = read_json(wd.store / "store.json")
    except Exception:
        return False
    if [s for s, _ in samples] != list(meta.get("samples", [])):
        return False
    S = meta.get("n_samples", 0)
    J = meta.get("n_junctions", 0)
    NS = meta.get("n_sites", 0)
    if not (S and J):
        return False
    from .store import STORE_FORMAT
    if int(meta.get("store_format", 0)) < STORE_FORMAT:
        log("Matrices construites par une version < 2.6 (comptages mal places "
            "dans la matrice) : 'build' refait, sans recompter")
        return False
    for name, ncol in (("K.split.dat", J), ("K.nonsplit.dat", NS)):
        p = winpath(wd.store / name)
        if not os.path.exists(p) or os.path.getsize(p) != S * ncol * 4:
            return False
    # 2.5 : un store construit avec le regroupement inexact (< 2.5) est refait
    try:
        from .store import check_blocks
        z = np.load(winpath(wd.store / "junctions.npz"))
        if check_blocks(z["jstart"], z["jend"], z["block_id"]):
            log("Matrices construites par une version < 2.5 (sites coupes entre "
                "blocs) : 'build' refait")
            return False
    except Exception:
        return False
    return True


def _analysis_ready(wd: Workdir, metrics) -> bool:
    """L'analyse est-elle deja faite pour toutes les metriques demandees ?"""
    try:
        meta = read_json(wd.store / "analysis.json")
    except Exception:
        return False
    S = meta.get("S", 0) or meta.get("n_samples", 0)
    if not S:
        return False
    try:   # matrices reconstruites depuis : l'analyse est perimee
        if os.path.getmtime(winpath(wd.store / "junctions.npz")) > \
                os.path.getmtime(winpath(wd.store / "analysis.json")):
            return False
    except OSError:
        return False
    for m in metrics:
        if m not in meta.get("info", {}):
            return False
        ncol = meta["info"][m]["n_cols"]
        for pre in ("pv_", "delta_"):
            p = winpath(wd.store / f"{pre}{m}.dat")
            if not os.path.exists(p) or os.path.getsize(p) != S * ncol * 4:
                return False
    return True


def cmd_all(a) -> None:
    t0 = time.time()
    missing = cmd_count(a)
    if missing:
        die("arret avant 'build' : des echantillons manquants laisseraient des "
            "lignes vides dans les matrices et fausseraient la cohorte.\n"
            "  - disque deconnecte : rebranchez-le et relancez la meme commande "
            "(le deja fait est garde) ;\n"
            "  - fichier abime ou d'un autre genome : retirez-le de "
            f"{Path(a.workdir) / 'samples.tsv'} puis relancez avec --samples.",
            code=1)
    wd = Workdir(a.workdir)
    samples = read_sample_sheet(str(wd.root / "samples.tsv"))
    metrics = tuple(a.metrics.split(","))
    if not a.force and _store_ready(wd, samples):
        log("Matrices deja construites : etape 'build' sautee (--force pour refaire)")
    else:
        cmd_build(a)
    # la reference doit etre fixee AVANT l'analyse ; si elle a change depuis
    # la derniere analyse, celle-ci est refaite
    refp = wd.store / "reference.json"
    before = refp.read_text(encoding="utf-8") if refp.exists() else ""
    _apply_reference(a, wd)
    after = refp.read_text(encoding="utf-8") if refp.exists() else ""
    if not a.force and _analysis_ready(wd, metrics) and before == after:
        log("Tests deja calcules : etape 'analyze' sautee (--force pour refaire)")
    else:
        cmd_analyze(a)
    if getattr(a, "retention", False):
        from .coverage import measure_coverage, save_coverage
        _, ann = _load_project(wd)
        meta = read_json(wd.store / "project.json")
        tg = _cov_targets(a, ann, meta["ref_names"], meta["ref_lengths"])
        S, M = measure_coverage(wd, samples, tg, meta["ref_lengths"],
                                threads=a.threads, chunk_mb=a.chunk_mb,
                                jobs=a.jobs, engine=a.engine,
                                ref_names=meta["ref_names"])
        save_coverage(wd, tg, samples, S, M)
    cmd_report(a)
    log(f"Pipeline complet termine en {(time.time()-t0)/60:.1f} min")


def cmd_bench(a) -> None:
    """Mesure la vitesse de lecture pour estimer la duree totale."""
    from .bam import scan_chunk
    from .bgzf import split_bam, vsplit

    samples = _samples_from_args(a)
    if not samples:
        die("indiquer au moins un BAM")
    sid, bam = samples[0]
    ref_names, ref_lengths = bam_refs(bam)
    ro = ref_offsets(ref_lengths)
    tot = sum(os.path.getsize(winpath(p)) for _, p in samples
              if os.path.exists(winpath(p)))
    from .scanner import resolve_engine
    if resolve_engine(a.engine) == "cpp":
        _bench_cpp(a, samples, bam, ro, tot)
        if getattr(a, "stress", 0):
            _stress_test(samples, ro, a)
        return
    sites = None
    if a.gtf:
        wd = Workdir(a.workdir).create()
        ann, *_ = _project(wd, samples, a.gtf)
        sites = ann.site_pos
    chunks = split_bam(bam, 64, min_chunk_bytes=a.bench_mb << 20)
    if len(chunks) < 2:
        die("BAM trop petit ou index .bai absent")
    v0, v1 = chunks[len(chunks) // 2]
    nbytes = (vsplit(v1)[0] - vsplit(v0)[0]) if v1 else 0
    t0 = time.time()
    res = scan_chunk(bam, v0, v1, ro, sites, anchor=a.anchor)
    dt = time.time() - t0
    speed = nbytes / max(dt, 1e-9)
    log(f"Tranche de {human(nbytes)} lue en {dt:.1f} s "
        f"({human(speed)}/s pour 1 coeur, {res.n_reads:,} reads)")
    log(f"Estimation : {human(tot)} au total -> "
        f"{tot/max(speed,1)/3600:.1f} h sur 1 coeur, "
        f"~{tot/max(speed,1)/3600/max(1,a.threads):.1f} h sur {a.threads} coeurs "
        f"(si le disque suit)")
    if getattr(a, "stress", 0):
        _stress_test(samples, ro, a)


def _bench_cpp(a, samples, bam, ro, tot) -> None:
    """Debit reel du lecteur C++ sur les 3 premiers millions de reads."""
    import tempfile
    from .scanner import run_sfscan, split_threads, write_ref_files

    with tempfile.TemporaryDirectory() as td:
        refs = write_ref_files(Path(td), [], ro)
        for t in sorted({1, max(1, a.threads)}):
            info = run_sfscan(["--bam", bam, "--ref-off", refs["ref_off"],
                               "-@", t, "--count", str(Path(td) / "b"),
                               "--max-reads", 3_000_000], retries=1,
                              label="bench")
            nb = info.get("compressed_bytes", 0)
            dt = max(info.get("seconds", 0.0), 1e-6)
            speed = nb / dt
            log(f"C++ -@ {t} : {info.get('records', 0):,} reads, {human(nb)} "
                f"en {dt:.1f} s -> {human(speed)}/s")
        jobs, per = split_threads(a.threads, a.jobs)
        log(f"Estimation : {human(tot)} de BAM -> ~{tot/max(speed,1)/3600:.1f} h "
            f"a -t {a.threads}, count ET retention en une seule passe "
            f"(`all --retention`), si le disque suit. Sur SSD, -j 2 a 4 "
            f"augmente encore le debit global.")


def _stress_test(samples, ro, a) -> None:
    """Lit en boucle plusieurs BAM pour voir si le disque tient la charge.

    Un disque externe vieillissant lache typiquement apres quelques minutes de
    lecture soutenue (alimentation, cable, controleur qui chauffe), pas des la
    premiere seconde : d'ou ce test de duree.
    """
    from concurrent.futures import ProcessPoolExecutor
    from .bam import scan_chunk
    from .bgzf import split_bam, vsplit
    from .count import worker
    from .util import io_incidents, reset_io_incidents

    minutes = float(a.stress)
    log(f"--- test d'endurance du disque : {minutes:.0f} min a {a.threads} "
        f"lecteur(s) ---")
    reset_io_incidents()
    tasks = []
    for sid, bam in samples:
        try:
            ch = split_bam(bam, 8, min_chunk_bytes=(a.bench_mb << 20))
        except Exception as e:  # noqa: BLE001
            warn(f"{sid} : index illisible ({e})")
            continue
        for v0, v1 in ch:
            tasks.append((bam, v0, v1, ro.tolist(), None,
                          dict(anchor=a.anchor)))
    if not tasks:
        die("aucune tranche lisible")
    t0 = time.time()
    nbytes = 0
    nchunk = 0
    last_log = -30.0
    ex = ProcessPoolExecutor(max_workers=max(1, a.threads)) if a.threads > 1 else None
    i = 0
    try:
        while time.time() - t0 < minutes * 60:
            batch = [tasks[(i + k) % len(tasks)] for k in range(max(1, a.threads))]
            i += len(batch)
            try:
                if ex is None:
                    outs = [worker(t) for t in batch]
                else:
                    outs = list(ex.map(worker, batch))
            except Exception as e:  # noqa: BLE001
                warn(f"ECHEC de lecture non rattrape : {e.__class__.__name__}: {e}")
                ex = ProcessPoolExecutor(max_workers=max(1, a.threads))
                continue
            for t, o in zip(batch, outs):
                c0 = vsplit(t[1])[0]
                c1 = vsplit(t[2])[0] if t[2] else c0
                nbytes += max(0, c1 - c0)
                nchunk += 1
            el = time.time() - t0
            if el - last_log < 30.0 and el < minutes * 60:
                continue
            last_log = el
            log(f"  {el/60:5.1f} min : {nchunk} tranches, {human(nbytes)} lus "
                f"({human(nbytes/max(el,1))}/s), "
                f"{io_incidents()} incident(s) rattrape(s)")
    finally:
        if ex is not None:
            ex.shutdown(wait=False, cancel_futures=True)
    n = io_incidents()
    if n == 0:
        log("Disque STABLE sur la duree du test : vous pouvez lancer l'analyse.")
    else:
        warn(f"{n} coupure(s) detectee(s) en {minutes:.0f} min. Elles sont "
             f"rattrapees automatiquement, mais elles ralentissent le run. "
             f"Essayez : -t 1, un autre cable/port USB (arriere du PC, pas de "
             f"hub), et desactivez la suspension selective USB.")


def cmd_selftest(a) -> None:
    from .selftest import run_selftest

    if not run_selftest(Path(a.workdir), threads=a.threads):
        die("auto-test en echec", code=1)


# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="splicefast",
        description="Detection d'evenements d'epissage aberrants (type FRASER), "
                    "optimisee pour de gros runs (lecteur BAM en C++).")
    p.add_argument("--version", action="version", version=f"splicefast {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp, need_wd=True):
        if need_wd:
            sp.add_argument("-w", "--workdir", required=True,
                            help="dossier de travail (sur un disque avec de la place)")
        sp.add_argument("-q", "--quiet", action="store_true")

    def inputs(sp):
        sp.add_argument("--run", nargs="*", default=[],
                        help="dossier(s) de run a scanner (recursif)")
        sp.add_argument("--samples", help="feuille TSV : sample<TAB>bam")
        sp.add_argument("--bam", nargs="*", default=[], help="BAM individuels")
        sp.add_argument("--pattern", default="*.bam")
        sp.add_argument("--no-recursive", action="store_true")

    sp = sub.add_parser("discover", help="lister les BAM d'un ou plusieurs runs")
    common(sp, need_wd=False)
    inputs(sp)
    sp.add_argument("--out", help="ecrire une feuille d'echantillons TSV")
    sp.set_defaults(func=cmd_discover, workdir=None)

    def engine_args(sp):
        sp.add_argument("--engine", choices=("auto", "cpp", "python"),
                        default="auto",
                        help="lecteur de BAM : cpp (sfscan, rapide), python "
                             "(sans compilation) ou auto (cpp si compile)")
        sp.add_argument("-j", "--jobs", type=int, default=1,
                        help="BAM lus simultanement (moteur cpp) ; les -t "
                             "coeurs sont repartis entre eux. 1 sur disque "
                             "externe, 2-4 sur SSD/NVMe")

    def retention_target_args(sp):
        sp.add_argument("--genes", nargs="*", default=[],
                        help="retention : symboles de genes a mesurer "
                             "(ex: SMC1A SETX ; defaut : tous)")
        sp.add_argument("--genes-file", help="fichier de symboles (un par ligne)")
        sp.add_argument("--edge", type=int, default=10,
                        help="retention : pb ignorees a chaque bord de l'intron (10)")
        sp.add_argument("--border-width", type=int, default=25,
                        help="retention : bases exoniques bordantes moyennees "
                             "comme reference (25 ; la base unique reste "
                             "mesuree a part)")
        sp.add_argument("--min-intron-len", type=int, default=80)

    def count_args(sp):
        inputs(sp)
        sp.add_argument("--gtf", help="annotation GTF/GTF.gz (fortement conseille)")
        sp.add_argument("--biotypes", default=None,
                        help="ne garder que ces gene_type (ex: protein_coding)")
        sp.add_argument("-t", "--threads", type=int,
                        default=max(1, (os.cpu_count() or 2) - 1))
        sp.add_argument("--anchor", type=int, default=5,
                        help="bases requises de chaque cote du site pour un "
                             "read non coupe (defaut 5, comme FRASER)")
        sp.add_argument("--min-anchor-split", type=int, default=5,
                        help="bases exoniques requises de chaque cote d'une jonction")
        sp.add_argument("--min-mapq", type=int, default=0)
        sp.add_argument("--min-intron", type=int, default=20)
        sp.add_argument("--max-intron", type=int, default=1000000)
        sp.add_argument("--keep-duplicates", action="store_true")
        sp.add_argument("--io-retries", type=int, default=20,
                        help="nombre de reprises sur erreur disque (disque "
                             "externe qui se deconnecte)")
        sp.add_argument("--io-wait", type=float, default=3.0,
                        help="attente initiale avant reprise, en secondes "
                             "(doublee a chaque echec)")
        sp.add_argument("--chunk-mb", type=int, default=32,
                        help="taille min d'une tranche de BAM par coeur (Mo)")
        sp.add_argument("--force", action="store_true",
                        help="recompter meme si deja fait")
        sp.add_argument("--force-annotation", action="store_true")
        engine_args(sp)
        sp.add_argument("--retention", action="store_true",
                        help="mesurer AUSSI la couverture des introns pendant "
                             "la meme lecture des BAM (module retention sans "
                             "seconde passe ; moteur C++ requis pour le gain)")
        retention_target_args(sp)

    sp = sub.add_parser("count", help="etape 1 : comptage par echantillon")
    common(sp)
    count_args(sp)
    sp.set_defaults(func=cmd_count)

    def reference_args(sp):
        g = sp.add_argument_group(
            "echantillons de reference",
            "Un echantillon hors reference est teste comme les autres mais ne "
            "participe pas a la norme (facteurs latents, valeur attendue, "
            "dispersion, medianes). Pour un sequencage trop peu profond qui "
            "sera repasse, ou tout echantillon dont on ne veut pas qu'il "
            "pollue la cohorte.")
        g.add_argument("--min-reads", type=float, default=0,
                       help="sous ce nombre de reads comptes, l'echantillon "
                            "est teste mais ne sert pas de reference "
                            "(ex: 60e6). 0 = pas de seuil")
        g.add_argument("--test-only", nargs="+", action="extend", default=[],
                       help="noms d'echantillons a tester sans les mettre "
                            "dans la reference")
        g.add_argument("--test-only-file",
                       help="fichier de noms (un par ligne) a exclure de la "
                            "reference")
        g.add_argument("--reference-match",
                       help="seuls les echantillons dont le NOM contient ce "
                            "motif (expression reguliere, sans casse) servent "
                            "de reference, ex: PAX ; les autres sont testes "
                            "seulement. Se combine avec --min-reads")
        g.add_argument("--reference-file",
                       help="fichier de noms : SEULS ceux-ci servent de "
                            "reference (tous les autres sont testes seuls)")

    def build_args(sp):
        sp.add_argument("--min-total", type=int, default=5,
                        help="reads totaux minimum sur la cohorte pour garder "
                             "une jonction")
        sp.add_argument("--min-samples", type=int, default=1)
        sp.add_argument("--block-size", type=int, default=20000)

    sp = sub.add_parser("build", help="etape 2 : construction des matrices")
    common(sp)
    build_args(sp)
    sp.add_argument("--samples", help="feuille TSV (defaut: <workdir>/samples.tsv)")
    reference_args(sp)
    sp.set_defaults(func=cmd_build)

    sp = sub.add_parser("panelapp",
                        help="telecharger PanelApp Australie (tous les panels, "
                             "avec la confiance de chaque gene)")
    common(sp)
    sp.add_argument("--update", action="store_true",
                    help="telecharger depuis l'API (defaut)")
    sp.add_argument("--from-json", help="importer des pages JSON deja telechargees")
    sp.add_argument("--save-json", help="garder aussi les pages JSON brutes ici")
    sp.add_argument("--out", help="fichier cache (defaut : <workdir>/panelapp/...)")
    sp.add_argument("--api", default="https://panelapp-aus.org/api/v1/genes/")
    sp.add_argument("--insecure", action="store_true",
                    help="ne pas verifier le certificat SSL (proxy hospitalier)")
    sp.set_defaults(func=cmd_panelapp)

    sp = sub.add_parser("probe",
                        help="relire une region dans le BAM et la comparer au "
                             "comptage stocke (evenement vu dans IGV mais absent)")
    common(sp)
    sp.add_argument("--sample", help="echantillon (le BAM est retrouve dans le projet)")
    sp.add_argument("--region", required=True, help="ex: chr5:37025000-37028000")
    sp.add_argument("--bam", help="ou un BAM explicite")
    sp.add_argument("--min-mapq", type=int, default=0)
    sp.add_argument("--anchor", type=int, default=5)
    sp.set_defaults(func=cmd_probe)

    sp = sub.add_parser("why",
                        help="pourquoi un evenement est-il absent du rapport ?")
    common(sp)
    sp.add_argument("--sample", required=True, help="echantillon a inspecter")
    sp.add_argument("--gene", help="symbole du gene (ex: DYNC1H1)")
    sp.add_argument("--chrom", help="ou une position : chromosome")
    sp.add_argument("--start", type=int, help="debut (1-based)")
    sp.add_argument("--end", type=int, help="fin")
    sp.add_argument("--padj", type=float, default=0.05)
    sp.add_argument("--delta", type=float, default=0.1)
    sp.add_argument("--nominal-p", type=float, default=1e-4)
    sp.add_argument("--top-per-sample", type=int, default=20)
    sp.add_argument("--top", type=int, default=300)
    sp.set_defaults(func=cmd_why)

    sp = sub.add_parser("reference",
                        help="choisir les echantillons qui definissent la "
                             "norme (sans recompter)")
    common(sp)
    reference_args(sp)
    sp.add_argument("--list", action="store_true",
                    help="lister les echantillons et leur statut")
    sp.set_defaults(func=cmd_reference)

    def analyze_args(sp):
        sp.add_argument("--metrics", default="jaccard,psi5,psi3,theta")
        sp.add_argument("--n-factors", default="auto",
                        help="nombre de facteurs latents (auto par defaut)")
        sp.add_argument("--sub-rows", type=int, default=30000)
        sp.add_argument("--min-k-one", type=int, default=10)
        sp.add_argument("--min-n-quantile", type=int, default=10)
        sp.add_argument("--min-cov", type=int, default=5)
        sp.add_argument("--robust-iter", type=int, default=3)
        sp.add_argument("--screen", type=float, default=0.1,
                        help="pre-filtre : test exact seulement si la p-value "
                             "approchee est sous ce seuil (0.1 par defaut ; "
                             "mettre 1.0 pour tout calculer exactement)")
        sp.add_argument("--covariates",
                        help="TSV sample<TAB>col1... (sexe, run, RIN...)")

    sp = sub.add_parser("analyze", help="etape 3 : modele beta-binomial")
    common(sp)
    analyze_args(sp)
    sp.set_defaults(func=cmd_analyze)

    def report_args(sp):
        sp.add_argument("--padj", type=float, default=0.05)
        sp.add_argument("--delta", type=float, default=0.1)
        sp.add_argument("--top", type=int, default=300,
                        help="max d'evenements par echantillon et par metrique")
        sp.add_argument("--sample", nargs="*", default=[],
                        help="limiter a certains echantillons")
        sp.add_argument("--shift-window", type=int, default=20)
        sp.add_argument("--pseudo-exon-max", type=int, default=1000)
        sp.add_argument("--panel", "--panels", nargs="+", action="extend",
                        default=[],
                        help="export(s) TSV PanelApp (Incidentalome.tsv, "
                             "Mendeliome.tsv...) : ajoute transmission et "
                             "phenotypes aux resultats")
        sp.add_argument("--html-rows", type=int, default=15000,
                        help=argparse.SUPPRESS)   # 2.5 : result.html remplace
        sp.add_argument("--panelapp", default=None,
                        help="cache PanelApp (defaut : <workdir>/panelapp/"
                             "panelapp_genes.tsv.gz, cree par `panelapp --update`)")
        sp.add_argument("--recurrent", type=int, default=4,
                        help="evenement 'recurrent' a partir de N echantillons "
                             "porteurs dans la cohorte (4)")
        sp.add_argument("--no-excel", action="store_true",
                        help="ne pas ecrire results.xlsx (gros runs)")
        sp.add_argument("--top-per-sample", type=int, default=20,
                        help="candidats de JONCTION a relire par echantillon "
                             "meme sans atteindre le FDR genome entier (20 ; "
                             "0 = FDR strict). Un evenement franc dans un gene "
                             "de panel passe souvent par la")
        sp.add_argument("--nominal-p", type=float, default=1e-4,
                        help="p brute maximale pour ces candidats (1e-4)")
        sp.add_argument("--no-harmonise", action="store_true",
                        help="ne pas regrouper les signaux en anomalies "
                             "(ni p combinee de retention)")
        sp.add_argument("--cov-top", type=int, default=20,
                        help="candidats de couverture a relire par "
                             "echantillon, meme sans atteindre padj (20 ; "
                             "0 pour s'en tenir au FDR strict)")
        sp.add_argument("--cov-min-samples", type=int, default=25,
                        help="echantillons couverts minimaux pour qu'un intron "
                             "soit testable (25)")
        sp.add_argument("--cov-min-border", type=float, default=10.0,
                        help="profondeur minimale aux bases bordantes pour "
                             "qu'un intron soit testable (module couverture)")
        sp.add_argument("--cov-min-ratio", type=float, default=0.1,
                        help="psi minimal (D/(D+bordures)) pour retenir un "
                             "evenement ; 0.1 ~ 11%% de retention, 0.33 ~ 50%%")

    sp = sub.add_parser("report", help="etape 4 : FDR + rapports")
    common(sp)
    report_args(sp)
    sp.set_defaults(func=cmd_report)

    sp = sub.add_parser("retention",
                        help="module dedie : retention d'intron par la couverture")
    common(sp)
    sp.add_argument("-t", "--threads", type=int, default=4)
    engine_args(sp)
    retention_target_args(sp)
    sp.add_argument("--panel", "--panels", nargs="+", action="extend",
                    default=[],
                    help="mesurer tous les genes de ces panels PanelApp")
    sp.add_argument("--all-genes", action="store_true",
                    help="(defaut) tous les introns annotes")
    sp.add_argument("--chunk-mb", type=int, default=64,
                    help="moteur python uniquement")
    sp.set_defaults(func=cmd_retention)

    sp = sub.add_parser("gene",
                        help="tableau couverture d'un gene, tous echantillons")
    common(sp)
    sp.add_argument("--gene", required=True, help="symbole du gene (ex: SMC1A)")
    sp.add_argument("--sample", help="echantillon a mettre en evidence")
    sp.add_argument("--top", type=int, default=8,
                    help="lignes affichees par intron (le TSV contient tout)")
    sp.add_argument("--min-border", type=float, default=10.0,
                    help="profondeur bordante minimale pour qu'un echantillon "
                         "soit exploitable (10)")
    sp.set_defaults(func=cmd_gene)

    sp = sub.add_parser("all", help="tout enchainer (count->build->analyze->report)")
    common(sp)
    count_args(sp)
    build_args(sp)
    reference_args(sp)
    analyze_args(sp)
    report_args(sp)
    sp.set_defaults(func=cmd_all)

    sp = sub.add_parser("bench", help="estimer la duree sur vos donnees")
    common(sp)
    inputs(sp)
    sp.add_argument("--gtf")
    sp.add_argument("--anchor", type=int, default=5)
    sp.add_argument("-t", "--threads", type=int,
                    default=max(1, (os.cpu_count() or 2) - 1))
    sp.add_argument("--bench-mb", type=int, default=64)
    engine_args(sp)
    sp.add_argument("--stress", type=float, default=0,
                    help="minutes de lecture continue pour tester la "
                         "stabilite du disque externe (ex: --stress 20)")
    sp.set_defaults(func=cmd_bench, biotypes=None, force_annotation=False)

    sp = sub.add_parser("selftest", help="test complet sur donnees simulees")
    common(sp)
    sp.add_argument("-t", "--threads", type=int, default=2)
    sp.set_defaults(func=cmd_selftest)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = build_parser()
    a = p.parse_args(argv)
    if getattr(a, "quiet", False):
        set_quiet(True)
    if getattr(a, "workdir", None):
        # journal persistant : indispensable pour un run de plusieurs heures
        try:
            from .util import set_logfile
            wdp = Path(a.workdir)
            (wdp / "logs").mkdir(parents=True, exist_ok=True)
            set_logfile(wdp / "logs" / "splicefast.log")
        except Exception:
            pass
    # duplication d'arguments entre sous-commandes (all) : valeurs par defaut
    for k, v in dict(min_total=5, min_samples=1, block_size=20000).items():
        if not hasattr(a, k):
            setattr(a, k, v)
    r = a.func(a)
    if a.cmd == "count" and r:
        return 1          # echantillons manquants : ne pas enchainer sur build
    return 0
