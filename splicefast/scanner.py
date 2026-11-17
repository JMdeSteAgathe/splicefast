"""
Pont vers le lecteur C++ `sfscan` (src/sfscan.cpp, compile par `make`).

Le binaire lit un BAM/CRAM entier en une passe, avec decompression
multi-thread (htslib), et produit au choix :
  * les comptages de l'etape `count` (jonctions + reads non coupes aux sites)
  * la couverture des introns du module `retention`
  * les deux a la fois, pour ne lire les BAM qu'une seule fois.

Les echanges se font par fichiers .npy : aucune extension Python a compiler,
et le binaire se teste seul en ligne de commande.

Moteurs :
  cpp    : sfscan obligatoire (erreur s'il n'est pas compile)
  python : l'ancien lecteur en Python pur (lent, mais sans compilation)
  auto   : cpp si le binaire est trouve, sinon python avec un avertissement
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from .util import die, log, wait_for_path, warn

PKG_BIN = Path(__file__).resolve().parent / "bin" / (
    "sfscan.exe" if os.name == "nt" else "sfscan")

# codes de sortie de sfscan
E_USAGE, E_HEADER, E_READ, E_WRITE, E_INPUT, E_CORRUPT = 2, 3, 4, 5, 6, 7
RETRYABLE = {E_READ, E_WRITE, -9, -15}   # disque, ecriture, processus tue

_WARNED = [False]


class ScanError(RuntimeError):
    def __init__(self, code: int, msg: str):
        super().__init__(msg)
        self.code = code


def find_sfscan() -> Optional[str]:
    env = os.environ.get("SPLICEFAST_SFSCAN")
    for cand in (env, str(PKG_BIN), shutil.which("sfscan")):
        if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return None


_EXE_OK: Dict[str, bool] = {}


def _exe_runs(exe: str) -> bool:
    """Le binaire se lance-t-il ? (antivirus, DLL manquante, autre OS...)"""
    if exe not in _EXE_OK:
        try:
            r = subprocess.run([exe, "--version"], stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True, timeout=60)
            _EXE_OK[exe] = r.returncode == 0 and "sfscan" in (r.stdout + r.stderr)
        except Exception as e:  # noqa: BLE001
            warn(f"sfscan ne peut pas etre lance ({e.__class__.__name__}: {e})")
            _EXE_OK[exe] = False
    return _EXE_OK[exe]


def resolve_engine(engine: str = "auto") -> str:
    engine = (engine or "auto").lower()
    if engine == "python":
        return "python"
    exe = find_sfscan()
    if exe and not _exe_runs(exe):
        if engine == "cpp":
            die(f"{exe} ne se lance pas (bloque par l'antivirus ?)")
        exe = None
    if exe:
        return "cpp"
    if engine == "cpp":
        die("lecteur C++ introuvable : compiler avec `make` a la racine du "
            "projet (voir README, section Installation), ou definir "
            "SPLICEFAST_SFSCAN=/chemin/vers/sfscan")
    if not _WARNED[0]:
        warn("lecteur C++ `sfscan` non compile : repli sur le lecteur Python "
             "(10 a 30 fois plus lent). Lancer `make` pour l'activer.")
        _WARNED[0] = True
    return "python"


def write_ref_files(directory: Path, ref_names: Sequence[str],
                    ref_off: np.ndarray) -> Dict[str, str]:
    """Noms + offsets des references du projet, verifies par sfscan sur
    chaque BAM (un BAM aligne sur un autre genome est refuse net au lieu de
    produire des comptages decales)."""
    directory.mkdir(parents=True, exist_ok=True)
    names = directory / "ref_names.txt"
    offs = directory / "ref_off.npy"
    # longueurs retrouvees depuis les offsets (pad de 1000, cf. bam.ref_offsets)
    ro = np.asarray(ref_off, dtype=np.int64)
    lens = (ro[1:] - ro[:-1] - 1000).tolist() if len(ro) > 1 else []
    txt = "".join(f"{n}\t{L}\n" for n, L in zip(ref_names, lens))
    if not names.exists() or names.read_text(encoding="utf-8") != txt:
        names.write_text(txt, encoding="utf-8")
    np.save(offs, np.asarray(ref_off, dtype=np.int64))
    return {"ref_names": str(names), "ref_off": str(offs)}


def _message(stderr: str, code: int) -> str:
    """Le message de sfscan (ligne 'sfscan: ...'), pas les avertissements htslib."""
    lines = [l.strip() for l in (stderr or "").splitlines() if l.strip()]
    own = [l for l in lines if l.startswith("sfscan:")]
    if own:
        return own[-1][len("sfscan:"):].strip()
    return lines[-1] if lines else f"code de sortie {code}"


def _file_state(path: Optional[str]):
    try:
        st = os.stat(path)
        with open(path, "rb") as fh:
            fh.read(1)
        return (st.st_size, st.st_mtime_ns)
    except (OSError, TypeError):
        return None


def run_sfscan(args: List[str], retries: int = 6, wait_path: Optional[str] = None,
               label: str = "") -> dict:
    """Lance sfscan ; reprend automatiquement sur erreur de lecture."""
    exe = find_sfscan()
    if exe is None:
        raise ScanError(E_USAGE, "sfscan introuvable")
    delay = 5.0
    state0 = _file_state(wait_path)
    for attempt in range(1, retries + 1):
        proc = subprocess.run([exe] + [str(a) for a in args],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True)
        if proc.returncode == 0:
            try:
                return json.loads(proc.stdout.strip().splitlines()[-1])
            except (ValueError, IndexError):
                return {}
        msg = _message(proc.stderr, proc.returncode)
        if proc.returncode not in RETRYABLE or attempt >= retries:
            raise ScanError(proc.returncode, msg)
        # Erreur de lecture alors que le fichier est lisible et inchange :
        # c'est le contenu qui est abime, pas le disque. Une seule relecture
        # immediate (secteur capricieux), puis on abandonne sans attendre.
        if proc.returncode == E_READ and state0 is not None \
                and _file_state(wait_path) == state0:
            if attempt >= 2:
                raise ScanError(E_CORRUPT, msg + " (fichier lisible mais "
                                "illisible deux fois de suite : corrompu)")
            warn(f"{label} : {msg} - relecture immediate")
            continue
        warn(f"{label} : {msg} - nouvelle tentative {attempt + 1}/{retries} "
             f"dans {delay:.0f}s")
        time.sleep(delay)
        if wait_path and not wait_for_path(wait_path, timeout=1800.0):
            raise ScanError(E_READ, f"disque inaccessible depuis 30 min : {wait_path}")
        delay = min(delay * 2, 120.0)
    raise ScanError(E_READ, "echec")


def split_threads(threads: int, jobs: int) -> tuple:
    """Repartit -t coeurs entre `jobs` BAM lus en parallele."""
    threads = max(1, int(threads))
    jobs = max(1, min(int(jobs or 1), threads))
    return jobs, max(1, threads // jobs)


def load_count_outputs(prefix: str, remove: bool = True) -> tuple:
    files = [f"{prefix}.jkey.npy", f"{prefix}.jcount.npy", f"{prefix}.nonsplit.npy"]
    keys, counts, ns = (np.load(f) for f in files)
    if remove:
        for f in files:
            try:
                os.remove(f)
            except OSError:
                pass
    return keys, counts, ns


def header_matches(bam: str, ref_names: Optional[Sequence[str]]) -> bool:
    """Entete identique au projet ? (exige par le lecteur Python, qui ne sait
    pas remapper les chromosomes ; le lecteur C++ le fait par nom)."""
    if not ref_names:
        return True
    from .bam import bam_refs
    names, _ = bam_refs(bam)
    return list(names) == list(ref_names)
