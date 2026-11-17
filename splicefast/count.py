"""
Etape 1 : comptage par echantillon.

Un echantillon = un BAM. On produit un .npz par echantillon :
    jkey   : cle entiere de la jonction (start global << 24 | longueur)
    jcount : nombre de reads coupes (K)
    nonsplit : reads NON coupes a chaque site d'epissage de la liste globale
    meta   : nb de reads, nb de reads epissés, duree

Points cles :
  * memoire bornee : on ne construit JAMAIS la matrice complete ici,
    seulement un dictionnaire de jonctions par echantillon (~10-50 Mo).
  * reprise sur incident : un echantillon deja compte est saute.
  * parallelisation intra-echantillon par tranches de BAM (index .bai),
    ce qui donne des lectures sequentielles par worker (bon pour un disque
    externe) et un usage CPU constant meme sur les gros BAM.
"""
from __future__ import annotations

import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import bam as bam_mod
from .bgzf import split_bam
from .util import (Progress, Workdir, human, io_incidents, log, read_json,
                   reset_io_incidents, retry_io, wait_for_path, warn, winpath,
                   write_json)

BAM_SUFFIX_STRIP = (
    ".Aligned.sortedByCoord.out", ".sortedByCoord.out", ".Aligned.out",
    ".sorted", ".dedup", ".markdup", ".recal", ".star", ".Aligned",
)

_SITES: Dict[str, np.ndarray] = {}


def _get_sites(path: Optional[str]) -> Optional[np.ndarray]:
    if not path:
        return None
    a = _SITES.get(path)
    if a is None:
        a = np.load(winpath(path))
        _SITES[path] = a
    return a


def worker(task):
    """Worker picklable (Windows -> spawn)."""
    bam_path, v0, v1, ref_off, sites_path, opts = task
    sites = _get_sites(sites_path)
    ro = np.asarray(ref_off, dtype=np.int64)
    res = bam_mod.scan_chunk(bam_path, v0, v1, ro, sites, **opts)
    keys = np.fromiter(res.junc.keys(), dtype=np.int64, count=len(res.junc))
    vals = np.fromiter(res.junc.values(), dtype=np.int32, count=len(res.junc))
    return keys, vals, res.nonsplit, res.n_reads, res.n_spliced


def clean_sample_id(p: Path) -> str:
    n = p.name
    for ext in (".bam", ".BAM"):
        if n.endswith(ext):
            n = n[: -len(ext)]
    for suf in BAM_SUFFIX_STRIP:
        if n.endswith(suf):
            n = n[: -len(suf)]
    return n


def discover_samples(run_dirs: Sequence[str], recursive: bool = True,
                     pattern: str = "*.bam") -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    seen: Dict[str, Path] = {}
    copies: List[Tuple[str, str]] = []
    for rd in run_dirs:
        root = Path(rd)
        if not root.exists():
            warn(f"Dossier introuvable : {rd}")
            continue
        it = root.rglob(pattern) if recursive else root.glob(pattern)
        # le plus proche de la racine d'abord : l'original avant ses copies
        for p in sorted(it, key=lambda x: (len(x.parts), str(x))):
            if p.name.endswith(".bai"):
                continue
            sid = clean_sample_id(p)
            if sid in seen:
                # Le meme BAM copie dans un sous-dossier (sauvegarde, tri par
                # tissu...) : le garder ferait de chaque patient son propre
                # "jumeau" dans la cohorte, et un evenement present chez deux
                # echantillons n'est plus aberrant -- il deviendrait invisible.
                if _same_file(seen[sid], p):
                    copies.append((str(p), str(seen[sid])))
                    continue
                sid = f"{p.parent.name}__{sid}"
                warn(f"deux BAM differents portent le meme nom : {seen.get(clean_sample_id(p))} "
                     f"et {p} -> le second est renomme {sid}")
            seen[sid] = p
            out.append((sid, str(p)))
    if copies:
        warn(f"{len(copies)} BAM ignores car copies identiques d'un BAM deja "
             f"retenu (ex. {copies[0][0]} = {copies[0][1]}). Un doublon dans la "
             f"cohorte masque les evenements du patient concerne.")
    return out


def _same_file(a: Path, b: Path) -> bool:
    """Meme contenu ? (taille + debut + fin : suffisant pour des BAM)."""
    try:
        if a.resolve() == b.resolve():
            return True
        sa, sb = a.stat().st_size, b.stat().st_size
        if sa != sb:
            return False
        with open(a, "rb") as fa, open(b, "rb") as fb:
            if fa.read(1 << 20) != fb.read(1 << 20):
                return False
            fa.seek(max(0, sa - (1 << 20)))
            fb.seek(max(0, sb - (1 << 20)))
            return fa.read() == fb.read()
    except OSError:
        return False


def read_sample_sheet(path: str) -> List[Tuple[str, str]]:
    out = []
    with open(winpath(path), "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            f = line.split("\t") if "\t" in line else line.split(",")
            if len(f) == 1:
                p = Path(f[0])
                out.append((clean_sample_id(p), str(p)))
            else:
                if f[0].lower() in ("sample", "sampleid", "sample_id"):
                    continue
                out.append((f[0].strip(), f[1].strip()))
    return out


def write_sample_sheet(path: Path, samples: Sequence[Tuple[str, str]],
                       groups: Optional[Dict[str, str]] = None) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(winpath(path), "w", encoding="utf-8") as fh:
        fh.write("sample\tbam\tgroup\n")
        for sid, p in samples:
            g = (groups or {}).get(sid, "")
            fh.write(f"{sid}\t{p}\t{g}\n")


# --------------------------------------------------------------------------- #

class Pool:
    """Pool de processus qui se reconstruit s'il casse.

    Un disque externe qui se deconnecte peut tuer un worker (et donc casser
    tout le pool : BrokenProcessPool). On le recree au lieu d'abandonner.
    """

    def __init__(self, threads: int):
        self.threads = threads
        self.ex = ProcessPoolExecutor(max_workers=threads) if threads > 1 else None

    def reset(self) -> None:
        if self.ex is not None:
            try:
                self.ex.shutdown(wait=False, cancel_futures=True)
            except Exception:
                pass
            warn("redemarrage des processus de lecture")
            self.ex = ProcessPoolExecutor(max_workers=self.threads)

    def shutdown(self) -> None:
        if self.ex is not None:
            try:
                self.ex.shutdown(wait=True)
            except Exception:
                pass


def count_sample(bam_path: str, ref_off: np.ndarray, sites_path: Optional[str],
                 opts: dict, threads: int, executor=None,
                 min_chunk_mb: int = 32, max_attempts: int = 6) -> dict:
    """Compte un echantillon en parallelisant par tranches de BAM.

    Chaque tranche est rejouable : si le disque se deconnecte pendant la
    lecture d'une tranche, seule cette tranche est relue (un worker ne
    fusionne rien tant qu'il n'a pas fini, donc aucun double comptage).
    """
    chunks = retry_io(split_bam, bam_path, max(1, threads),
                      min_chunk_bytes=min_chunk_mb << 20,
                      what=f"decoupage de {os.path.basename(bam_path)}",
                      path=winpath(bam_path))
    tasks = [(bam_path, a, b, ref_off.tolist(), sites_path, opts) for a, b in chunks]
    key_parts: List[np.ndarray] = []
    val_parts: List[np.ndarray] = []
    nonsplit = None
    n_reads = n_spl = 0
    t0 = time.time()

    def merge(r):
        nonlocal nonsplit, n_reads, n_spl
        keys, vals, ns, nr, nsp = r
        if len(keys):
            key_parts.append(keys)
            val_parts.append(vals.astype(np.int64))
        if ns is not None:
            if nonsplit is None:
                nonsplit = ns.astype(np.int64)
            else:
                nonsplit += ns
        n_reads += nr
        n_spl += nsp

    pool = executor
    ex = getattr(pool, "ex", pool)
    remaining = tasks
    for attempt in range(1, max_attempts + 1):
        failed: List[tuple] = []
        last_err: Optional[BaseException] = None
        if ex is None or len(remaining) == 1:
            for t in remaining:
                try:
                    merge(worker(t))
                except Exception as e:  # noqa: BLE001
                    failed.append(t)
                    last_err = e
        else:
            futs = {ex.submit(worker, t): t for t in remaining}
            broken = False
            for f in as_completed(futs):
                try:
                    merge(f.result())
                except BrokenProcessPool as e:
                    failed.append(futs[f])
                    last_err = e
                    broken = True
                except Exception as e:  # noqa: BLE001
                    failed.append(futs[f])
                    last_err = e
            if broken and pool is not None and hasattr(pool, "reset"):
                pool.reset()
                ex = pool.ex
        if not failed:
            break
        if attempt >= max_attempts:
            raise RuntimeError(
                f"{os.path.basename(bam_path)} : {len(failed)} tranche(s) "
                f"illisibles apres {max_attempts} tentatives") from last_err
        warn(f"{os.path.basename(bam_path)} : {len(failed)} tranche(s) a relire "
             f"(tentative {attempt + 1}/{max_attempts}) - derniere erreur : "
             f"{last_err.__class__.__name__}: {last_err}")
        if not wait_for_path(bam_path, timeout=1800.0):
            raise RuntimeError(
                f"disque inaccessible depuis 30 min : {bam_path}") from last_err
        # on repart de zero pour cet echantillon si des donnees ont pu etre
        # fusionnees partiellement lors d'une tentative precedente
        remaining = failed

    if key_parts:
        ak = np.concatenate(key_parts)
        av = np.concatenate(val_parts)
        del key_parts[:], val_parts[:]
        uk, inv = np.unique(ak, return_inverse=True)
        sums = np.bincount(inv, weights=av, minlength=len(uk)).astype(np.int64)
    else:
        uk = np.zeros(0, np.int64)
        sums = np.zeros(0, np.int64)
    dt = time.time() - t0
    return dict(jkey=uk, jcount=sums, nonsplit=nonsplit, n_reads=n_reads,
                n_spliced=n_spl, seconds=dt, n_chunks=len(chunks))


def bam_signature(bam: str) -> str:
    """Identite d'un BAM : chemin + taille + date de modification.

    2.5 : stockee dans chaque fichier de comptage. Un BAM remplace (realigne,
    recopie, autre fichier portant le meme nom d'echantillon) est detecte et
    recompte au lieu de reutiliser silencieusement un comptage perime.
    """
    try:
        st = os.stat(winpath(bam))
        return f"{os.path.abspath(bam)}|{st.st_size}|{int(st.st_mtime)}"
    except OSError:
        return f"{bam}|?|?"


def _same_bam(old: str, new: str) -> bool:
    """Meme BAM ? La taille et la date suffisent : le chemin peut changer
    (disque monte ailleurs, lettre de lecteur differente sous Windows)."""
    a, b = old.split("|"), new.split("|")
    if len(a) < 3 or len(b) < 3 or "?" in a[1:] or "?" in b[1:]:
        return True
    return a[1] == b[1] and abs(int(a[2]) - int(b[2])) <= 2


def save_sample(path: Path, sid: str, res: dict, sites_hash: str, opts: dict,
                bam: str = "") -> None:
    keys = np.asarray(res["jkey"], dtype=np.int64)
    vals = np.minimum(np.asarray(res["jcount"]),
                      np.iinfo(np.int32).max).astype(np.int32)
    ns = res["nonsplit"]
    ns = (np.zeros(0, np.int32) if ns is None
          else np.minimum(ns, np.iinfo(np.int32).max).astype(np.int32))
    tmp = path.with_suffix(".part")

    def _write():
        with open(winpath(tmp), "wb") as fh:
            np.savez(
                fh, jkey=keys, jcount=vals, nonsplit=ns,
                meta=np.array([res["n_reads"], res["n_spliced"], len(keys),
                               int(res["seconds"]), res["n_chunks"]],
                              dtype=np.int64),
                info=np.array([sid, sites_hash, repr(opts),
                               bam_signature(bam) if bam else ""], dtype=object),
            )
        os.replace(winpath(tmp), winpath(path))

    retry_io(_write, what=f"ecriture du comptage de {sid}",
             path=str(path.parent))


def count_all(wd: Workdir, samples: Sequence[Tuple[str, str]], ref_off: np.ndarray,
              sites_path: Optional[str], sites_hash: str, opts: dict,
              threads: int = 4, force: bool = False,
              min_chunk_mb: int = 32, engine: str = "auto", jobs: int = 1,
              ref_names: Optional[Sequence[str]] = None,
              cov_targets=None) -> List[str]:
    """Comptage de tous les echantillons (reprenable).

    engine="cpp" : lecteur C++ sfscan (voir scanner.py). Si `cov_targets` est
    fourni, la couverture des introns est mesuree PENDANT la meme lecture et
    mise en cache : le module `retention` n'aura plus a relire les BAM.
    """
    from .scanner import header_matches, resolve_engine
    if resolve_engine(engine) == "cpp":
        return count_all_cpp(wd, samples, ref_off, sites_path, sites_hash, opts,
                             threads=threads, force=force, jobs=jobs,
                             ref_names=ref_names, cov_targets=cov_targets)
    bad = [sid for sid, p in samples if not header_matches(p, ref_names)]
    if bad:
        warn(f"{len(bad)} BAM avec des chromosomes differents du projet "
             f"({', '.join(bad[:5])}...) : le lecteur Python ne sait pas les "
             f"remapper, ils sont ecartes. Compiler sfscan (make) pour les "
             f"traiter.")
        samples = [x for x in samples if x[0] not in set(bad)]
    todo = []
    for sid, p in samples:
        out = wd.sample_counts(sid)
        if out.exists() and not force and _cache_valid(out, sid, p, sites_hash):
            continue
        todo.append((sid, p))
    log(f"Comptage : {len(todo)} echantillon(s) a traiter "
        f"({len(samples) - len(todo)} deja fait(s))")
    if not todo:
        return []

    total_bytes = 0
    for _, p in todo:
        try:
            total_bytes += os.path.getsize(winpath(p))
        except OSError:
            pass
    log(f"Volume BAM a lire : {human(total_bytes)}")

    done_bytes = 0
    t_start = time.time()
    reset_io_incidents()
    pool = Pool(threads)
    skipped: List[str] = []
    try:
        for i, (sid, p) in enumerate(todo, 1):
            try:
                sz = os.path.getsize(winpath(p)) if os.path.exists(winpath(p)) else 0
            except OSError:
                sz = 0
            before = io_incidents()
            try:
                res = count_sample(p, ref_off, sites_path, opts, threads, pool,
                                   min_chunk_mb=min_chunk_mb)
            except Exception as e:  # noqa: BLE001
                # on n'abandonne pas tout le run pour un echantillon : on le
                # met de cote et on continue. Relancer la commande le reprendra.
                warn(f"{sid} : ECHEC ({e.__class__.__name__}: {e}) - "
                     f"echantillon mis de cote, on continue")
                skipped.append(sid)
                pool.reset()
                if not wait_for_path(p, timeout=1800.0):
                    warn("le disque ne repond plus : arret propre. "
                         "Rebranchez-le et relancez la meme commande.")
                    break
                continue
            save_sample(wd.sample_counts(sid), sid, res, sites_hash, opts, bam=p)
            done_bytes += sz
            inc = io_incidents() - before
            el = time.time() - t_start
            speed = done_bytes / max(el, 1e-9)
            eta = (total_bytes - done_bytes) / max(speed, 1e-9)
            log(f"[{i}/{len(todo)}] {sid} : {res['n_reads']:,} reads, "
                f"{len(res['jkey']):,} jonctions, {res['seconds']/60:.1f} min "
                f"({human(sz)}, {human(speed)}/s) - reste ~{eta/3600:.1f} h"
                + (f" [{inc} incident(s) disque rattrape(s)]" if inc else ""))
    finally:
        pool.shutdown()
    if io_incidents():
        warn(f"{io_incidents()} incident(s) disque rattrape(s) pendant le "
             f"comptage. Si cela se repete, baissez -t (essayez -t 2) et "
             f"verifiez l'alimentation/le cable du disque externe.")
    if skipped:
        warn(f"{len(skipped)} echantillon(s) non comptes : "
             f"{', '.join(skipped[:8])}{' ...' if len(skipped) > 8 else ''}. "
             f"Relancez la meme commande pour les reprendre.")
    return skipped


# --------------------------------------------------------------------------- #
# Moteur C++ (sfscan)
# --------------------------------------------------------------------------- #
def _todo(wd: Workdir, samples, sites_hash: str, force: bool):
    todo = []
    for sid, p in samples:
        out = wd.sample_counts(sid)
        if out.exists() and not force and _cache_valid(out, sid, p, sites_hash):
            continue
        todo.append((sid, p))
    return todo


def _cache_valid(out: Path, sid: str, bam: str, sites_hash: str) -> bool:
    try:
        z = np.load(winpath(out), allow_pickle=True)
        info = z["info"]
    except Exception:
        warn(f"{sid} : fichier de comptage illisible -> recomptage")
        return False
    if str(info[1]) != sites_hash:
        warn(f"{sid} : liste de sites differente -> recomptage")
        return False
    if len(info) >= 4 and str(info[3]):
        if not _same_bam(str(info[3]), bam_signature(bam)):
            warn(f"{sid} : le BAM a change depuis le comptage "
                 f"({str(info[3]).split('|')[0]} -> {bam}) -> recomptage")
            return False
    else:
        log(f"{sid} : comptage anterieur a la 2.5 (BAM d'origine non "
            f"enregistre) : reutilise. --force pour recompter.")
    return True


def count_all_cpp(wd: Workdir, samples: Sequence[Tuple[str, str]],
                  ref_off: np.ndarray, sites_path: Optional[str],
                  sites_hash: str, opts: dict, threads: int = 4,
                  force: bool = False, jobs: int = 1,
                  ref_names: Optional[Sequence[str]] = None,
                  cov_targets=None) -> List[str]:
    from concurrent.futures import ThreadPoolExecutor
    from .bam import DEFAULT_EXCLUDE
    from .coverage import (cov_cache_ok, cov_cache_save, cov_sample_hash,
                           cov_scan_args)
    from .scanner import (ScanError, load_count_outputs, run_sfscan,
                          split_threads, write_ref_files)

    todo = _todo(wd, samples, sites_hash, force)
    log(f"Comptage (C++) : {len(todo)} echantillon(s) a traiter "
        f"({len(samples) - len(todo)} deja fait(s))")
    if not todo:
        return []
    jobs, per = split_threads(threads, jobs)
    total_bytes = sum(os.path.getsize(p) for _, p in todo if os.path.exists(p))
    log(f"Volume BAM a lire : {human(total_bytes)} - {jobs} BAM a la fois, "
        f"{per} thread(s) de decompression chacun")
    refs = write_ref_files(wd.tmp, ref_names or [], ref_off)
    scan_dir = wd.tmp / "scan"
    scan_dir.mkdir(parents=True, exist_ok=True)
    cov_mapq, cov_excl = 0, DEFAULT_EXCLUDE    # reglages du module retention

    def one(item):
        sid, bam = item
        from .util import safe_name
        prefix = str(scan_dir / safe_name(sid))
        args = ["--bam", bam, "--ref-off", refs["ref_off"], "-@", per,
                "--count", prefix,
                "--anchor", opts.get("anchor", 5),
                "--min-mapq", opts.get("min_mapq", 0),
                "--exclude", opts.get("exclude_flags", DEFAULT_EXCLUDE),
                "--min-intron", opts.get("min_intron", 20),
                "--max-intron", opts.get("max_intron", 1_000_000),
                "--min-anchor-split", opts.get("min_anchor_split", 5)]
        if ref_names:
            args += ["--ref-names", refs["ref_names"]]
        if sites_path:
            args += ["--sites", sites_path]
        cov_out = None
        if cov_targets is not None:
            ch = cov_sample_hash(cov_targets, bam, cov_mapq, cov_excl)
            if not cov_cache_ok(wd, sid, ch, cov_targets.n):
                extra, cov_out = cov_scan_args(wd, cov_targets, sid, cov_mapq,
                                               cov_excl)
                args += extra
        t0 = time.time()
        info = run_sfscan(args, wait_path=bam, label=sid)
        keys, counts, ns = load_count_outputs(prefix)
        res = dict(jkey=keys, jcount=counts,
                   nonsplit=(ns if sites_path else None),
                   n_reads=int(info.get("n_reads", 0)),
                   n_spliced=int(info.get("n_spliced", 0)),
                   seconds=time.time() - t0, n_chunks=1,
                   other_contig=int(info.get("other_contig_records", 0)))
        save_sample(wd.sample_counts(sid), sid, res, sites_hash, opts, bam=bam)
        if cov_out is not None:
            cov_cache_save(wd, sid, ch, np.load(cov_out))
            os.remove(cov_out)
        return res, cov_out is not None

    sizes = {sid: (os.path.getsize(p) if os.path.exists(p) else 0)
             for sid, p in todo}
    done_bytes = 0
    t_start = time.time()
    skipped: List[str] = []
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(one, it): it[0] for it in todo}
        for i, f in enumerate(as_completed(futs), 1):
            sid = futs[f]
            try:
                res, with_cov = f.result()
            except ScanError as e:
                warn(f"{sid} : ECHEC ({e}) - echantillon mis de cote, on continue")
                skipped.append(sid)
                continue
            done_bytes += sizes[sid]
            el = time.time() - t_start
            speed = done_bytes / max(el, 1e-9)
            eta = (total_bytes - done_bytes) / max(speed, 1e-9)
            other = res.get("other_contig", 0)
            log(f"[{i}/{len(todo)}] {sid} : {res['n_reads']:,} reads, "
                f"{len(res['jkey']):,} jonctions"
                + (" + couverture" if with_cov else "")
                + (f" ({other:,} reads sur des sequences hors projet, ignores)"
                   if other else "")
                + f", {res['seconds']/60:.1f} min ({human(sizes[sid])}, "
                f"{human(speed)}/s) - reste ~{eta/3600:.1f} h")
    if skipped:
        warn(f"{len(skipped)} echantillon(s) non comptes : "
             f"{', '.join(skipped[:8])}{' ...' if len(skipped) > 8 else ''}. "
             f"Relancez la meme commande pour les reprendre.")
    return skipped
