"""Utilitaires communs : logs, chemins, ecriture atomique, layout du workdir."""
from __future__ import annotations

import json
import os
import re
import sys
import time
from pathlib import Path
from typing import Optional

_T0 = time.time()
_QUIET = False


def set_quiet(q: bool) -> None:
    global _QUIET
    _QUIET = q


_LOGFILE = [None]


def set_logfile(path) -> None:
    """Duplique les messages dans un fichier journal (utile pour un run de
    plusieurs heures : on peut relire apres coup tous les incidents disque)."""
    try:
        fh = open(winpath(path), "a", encoding="utf-8", errors="replace")
        _LOGFILE[0] = fh
    except OSError:
        _LOGFILE[0] = None


def log(msg: str, *, level: str = "INFO") -> None:
    if _QUIET and level == "INFO":
        if _LOGFILE[0] is None:
            return
    el = time.time() - _T0
    h, rem = divmod(int(el), 3600)
    m, s = divmod(rem, 60)
    line = f"[{h:02d}:{m:02d}:{s:02d}] {level:<5s} {msg}\n"
    sys.stderr.write(line)
    sys.stderr.flush()
    fh = _LOGFILE[0]
    if fh is not None:
        try:  # best effort : le journal ne doit jamais faire planter le run
            fh.write(time.strftime("%Y-%m-%d %H:%M:%S ") + line)
            fh.flush()
        except Exception:
            pass


def warn(msg: str) -> None:
    log(msg, level="WARN")


def die(msg: str, code: int = 2) -> None:
    log(msg, level="ERROR")
    sys.exit(code)


def human(n: float) -> str:
    for unit in ("o", "Ko", "Mo", "Go", "To"):
        if abs(n) < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024.0
    return f"{n:.1f} Po"


def winpath(p: str | os.PathLike) -> str:
    """Chemin long-path-safe sous Windows (evite la limite MAX_PATH de 260)."""
    s = str(p)
    if os.name == "nt":
        s = os.path.abspath(s)
        if not s.startswith("\\\\?\\") and len(s) > 240:
            if s.startswith("\\\\"):
                s = "\\\\?\\UNC\\" + s[2:]
            else:
                s = "\\\\?\\" + s
    return s


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._+-]+", "_", s)


class Workdir:
    """Arborescence de travail. Tout ce qui est volumineux vit ici (ex: /data/splicefast)."""

    def __init__(self, root: str | os.PathLike):
        self.root = Path(root)
        self.counts = self.root / "counts"        # 1 fichier .npz par echantillon
        self.store = self.root / "store"          # matrices memmap + index
        self.results = self.root / "results"      # TSV / HTML
        self.tmp = self.root / "tmp"
        self.covcounts = self.root / "covcounts"  # couverture par echantillon
        self.logs = self.root / "logs"

    def create(self) -> "Workdir":
        for d in (self.root, self.counts, self.store, self.results, self.tmp,
                  self.logs, self.covcounts):
            d.mkdir(parents=True, exist_ok=True)
        # rediriger les temporaires Python sur le meme volume (pas de C:)
        os.environ["TMPDIR"] = str(self.tmp)
        os.environ["TEMP"] = str(self.tmp)
        os.environ["TMP"] = str(self.tmp)
        return self

    def sample_counts(self, sid: str) -> Path:
        return self.counts / f"{safe_name(sid)}.npz"

    def __repr__(self) -> str:
        return f"Workdir({self.root})"


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    with open(winpath(tmp), "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=1, ensure_ascii=False)
    os.replace(winpath(tmp), winpath(path))


def read_json(path: Path):
    with open(winpath(path), "r", encoding="utf-8") as fh:
        return json.load(fh)


def free_space(path: str | os.PathLike) -> int:
    try:
        import shutil

        return shutil.disk_usage(str(path)).free
    except Exception:
        return -1


def check_space(path: str | os.PathLike, needed: int) -> None:
    free = free_space(path)
    if 0 <= free < needed:
        warn(f"Espace libre faible sur {path} : {human(free)} (recommande >= {human(needed)})")


def open_maybe_gzip(path: str | os.PathLike, mode: str = "rt"):
    p = str(path)
    if p.endswith(".gz") or p.endswith(".bgz"):
        import gzip

        return gzip.open(winpath(p), mode, encoding="utf-8" if "t" in mode else None)
    return open(winpath(p), mode, encoding="utf-8" if "t" in mode else None)


class Progress:
    """Barre de progression minimaliste, sur stderr, sans dependance."""

    def __init__(self, total: int, label: str = "", every: float = 5.0):
        self.total = max(1, total)
        self.n = 0
        self.label = label
        self.every = every
        self.t0 = time.time()
        self.last = 0.0

    def update(self, k: int = 1) -> None:
        self.n += k
        now = time.time()
        if now - self.last < self.every and self.n < self.total:
            return
        self.last = now
        frac = self.n / self.total
        el = now - self.t0
        eta = el / max(frac, 1e-9) - el
        bar = "#" * int(frac * 28)
        sys.stderr.write(
            f"\r  {self.label} [{bar:<28s}] {100*frac:5.1f}%  {self.n}/{self.total}"
            f"  ecoule {el/60:.1f}min  reste ~{eta/60:.1f}min   "
        )
        sys.stderr.flush()

    def close(self) -> None:
        sys.stderr.write("\n")
        sys.stderr.flush()


# --------------------------------------------------------------------------- #
# Resistance aux disques externes qui se deconnectent
# --------------------------------------------------------------------------- #

IO_RETRIES = int(os.environ.get("SPLICEFAST_IO_RETRIES", "20"))
IO_WAIT = float(os.environ.get("SPLICEFAST_IO_WAIT", "3.0"))
IO_WAIT_MAX = float(os.environ.get("SPLICEFAST_IO_WAIT_MAX", "120.0"))

# compteur d'incidents, remis a zero par le processus principal
_IO_INCIDENTS = [0]


def io_incidents() -> int:
    return _IO_INCIDENTS[0]


def reset_io_incidents() -> None:
    _IO_INCIDENTS[0] = 0


def _is_transient(exc: BaseException) -> bool:
    """Une erreur disque qui a des chances de disparaitre si on reessaie.

    Sous Windows, un disque USB/eSATA qui se re-enumere fait remonter
    ERROR_INVALID_PARAMETER (87) -> EINVAL (errno 22) sur un handle deja
    ouvert, ainsi que ERROR_DEVICE_NOT_CONNECTED / ERROR_NOT_READY.
    """
    if isinstance(exc, (TimeoutError, InterruptedError)):
        return True
    if not isinstance(exc, OSError):
        return False
    import errno
    if exc.errno in (errno.EINVAL, errno.EIO, errno.ENODEV, errno.ENXIO,
                     errno.EBUSY, errno.EAGAIN, errno.ENOENT, errno.EACCES,
                     errno.EBADF, errno.ETIMEDOUT):
        return True
    return getattr(exc, "winerror", None) in (
        1, 21, 87, 433, 1117, 1167, 1392, 6, 995, 1006)


def wait_for_path(path: str | os.PathLike, timeout: float = 900.0,
                  poll: float = 5.0) -> bool:
    """Attend qu'un chemin redevienne lisible (disque rebranche / reveille)."""
    p = winpath(path)
    t0 = time.time()
    while time.time() - t0 < timeout:
        try:
            with open(p, "rb") as fh:
                fh.read(1)
            return True
        except OSError:
            time.sleep(poll)
    return False


def retry_io(fn, *args, what: str = "operation disque",
             retries: Optional[int] = None, wait: Optional[float] = None,
             path: Optional[str] = None, **kw):
    """Execute fn en reessayant sur les erreurs disque transitoires."""
    retries = IO_RETRIES if retries is None else retries
    delay = IO_WAIT if wait is None else wait
    for attempt in range(1, retries + 1):
        try:
            return fn(*args, **kw)
        except Exception as e:  # noqa: BLE001
            if attempt >= retries or not _is_transient(e):
                raise
            _IO_INCIDENTS[0] += 1
            warn(f"{what} : erreur disque ({e.__class__.__name__}: {e}) - "
                 f"nouvelle tentative {attempt}/{retries} dans {delay:.0f}s")
            time.sleep(delay)
            if path:
                wait_for_path(path, timeout=max(delay * 4, 60.0))
            delay = min(delay * 2, IO_WAIT_MAX)


class RobustReader:
    """Fichier binaire en lecture seule qui survit aux deconnexions du disque.

    On memorise nous-memes la position logique : quand le peripherique se
    re-enumere, le handle ouvert est definitivement mort (Windows renvoie
    EINVAL / errno 22), il faut donc rouvrir le fichier par son chemin et se
    repositionner. Vu d'en haut, la lecture reprend exactement ou elle en
    etait.
    """

    def __init__(self, path, buffering: int = 4 << 20,
                 retries: Optional[int] = None, wait: Optional[float] = None):
        self.path = winpath(path)
        self.buffering = buffering
        # lus a l'execution : les reglages peuvent avoir ete changes en ligne
        # de commande apres l'import du module
        self.retries = IO_RETRIES if retries is None else retries
        self.wait = IO_WAIT if wait is None else wait
        self.pos = 0
        self.fh = None
        self.incidents = 0
        self._open()

    def _open(self) -> None:
        self.fh = open(self.path, "rb", buffering=self.buffering)
        if self.pos:
            self.fh.seek(self.pos)

    def _drop(self) -> None:
        try:
            if self.fh is not None:
                self.fh.close()
        except Exception:
            pass
        self.fh = None

    def _recover(self, exc: BaseException, attempt: int) -> None:
        self.incidents += 1
        _IO_INCIDENTS[0] += 1
        self._drop()
        delay = min(self.wait * (2 ** (attempt - 1)), IO_WAIT_MAX)
        warn(f"lecture de {os.path.basename(self.path)} interrompue "
             f"({exc.__class__.__name__}: {exc}) - reprise a l'octet "
             f"{self.pos} dans {delay:.0f}s (essai {attempt}/{self.retries})")
        time.sleep(delay)
        wait_for_path(self.path, timeout=max(delay * 5, 120.0))

    def read(self, n: int) -> bytes:
        for attempt in range(1, self.retries + 1):
            try:
                if self.fh is None:
                    self._open()
                b = self.fh.read(n)
                self.pos += len(b)
                return b
            except Exception as e:  # noqa: BLE001
                if attempt >= self.retries or not _is_transient(e):
                    raise
                self._recover(e, attempt)
        raise OSError(f"{self.path}: disque inaccessible apres "
                      f"{self.retries} tentatives")

    def seek(self, offset: int, whence: int = 0) -> int:
        if whence != 0:
            raise ValueError("RobustReader : seek absolu uniquement")
        for attempt in range(1, self.retries + 1):
            try:
                if self.fh is None:
                    self._open()
                self.fh.seek(offset)
                self.pos = offset
                return offset
            except Exception as e:  # noqa: BLE001
                if attempt >= self.retries or not _is_transient(e):
                    raise
                self.pos = offset
                self._recover(e, attempt)
        raise OSError(f"{self.path}: disque inaccessible")

    def tell(self) -> int:
        return self.pos

    def close(self) -> None:
        self._drop()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()
