#!/usr/bin/env python3
"""autodE single-shot installer — driver.

Hardware-aware, resumable, rich-powered. Invoked by ``install.sh`` (which
bootstraps uv + micromamba and runs this under ``uv run --with rich``).

Design goals (see docs/handoff/autode-install-manifest.md):
  * ONE command sets up autodE correctly for the detected hardware.
  * Every step is idempotent and CHECKPOINTED — a failure never restarts from
    zero (``--resume`` skips completed steps; ``--fix STEP`` re-runs one).
  * Modern tooling: uv (pip/venv), micromamba (envs). SOTA logging via rich.
  * Optional method backends (UMA/GPU4PySCF/VeloxChem/g-xTB/MetalloGen) are
    detected + wired only where the hardware supports them.
  * External tools we can't install (ORCA) are PATH/version-checked with clear
    manual-install guidance; secrets (HF token for off-cluster UMA) are prompted.

Run ``./install.sh --help`` for options.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Callable, Optional

try:
    from rich.console import Console, Group
    from rich.panel import Panel
    from rich.table import Table
    from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn
    from rich.prompt import Confirm, Prompt
    from rich.text import Text
    from rich.rule import Rule
    console = Console()
    _RICH = True
except ImportError:
    # Graceful degradation: on a node with no egress (can't fetch rich) the
    # installer still runs with plain output. Minimal shims mirror the tiny
    # slice of the rich API this driver uses.
    _RICH = False

    class _Console:
        def print(self, *a, **k):
            for x in a:
                s = getattr(x, "_plain", None)
                print(s if s is not None else _strip(x))
        def rule(self, title=""):
            print(f"\n=== {_strip(title)} ===")
    console = _Console()

    def _strip(x):
        import re as _re
        return _re.sub(r"\[/?[a-z0-9 #._]*\]", "", str(x), flags=_re.I)

    class _Plain:
        def __init__(self, text="", *a, **k): self._plain = _strip(text)
        def __str__(self): return self._plain
        @staticmethod
        def assemble(*parts, **k):
            out = []
            for p in parts:
                out.append(_strip(p[0]) if isinstance(p, tuple) else _strip(p))
            return _Plain("".join(out))
    class Panel(_Plain):
        def __init__(self, body, title="", **k):
            super().__init__(f"[{_strip(title)}]\n{_strip(body)}" if title else _strip(body))
        @staticmethod
        def fit(body, title="", **k): return Panel(body, title)
    Text = _Plain
    def Group(*a, **k): return "\n".join(_strip(x) for x in a)
    def Rule(t=""): return f"--- {_strip(t)} ---"

    class Table:
        def __init__(self, *a, **k): self.rows = []
        def add_column(self, *a, **k): pass
        def add_row(self, *cells, **k): self.rows.append(" | ".join(_strip(c) for c in cells))
        @property
        def _plain(self): return "\n".join(self.rows)
    class Confirm:
        @staticmethod
        def ask(q, default=True):
            try: return (input(f"{_strip(q)} [{'Y/n' if default else 'y/N'}] ").strip().lower() or ("y" if default else "n")).startswith("y")
            except EOFError: return default
    class Prompt:
        @staticmethod
        def ask(q, default="", password=False):
            try: return input(f"{_strip(q)}: ").strip() or default
            except EOFError: return default
    class Progress:
        def __init__(self, *a, **k): pass
        def __enter__(self): return self
        def __exit__(self, *a): pass
        def add_task(self, desc, total=None): print(_strip(desc)); return 0
        def update(self, *a, **k): pass
        def advance(self, *a, **k): pass
    def SpinnerColumn(*a, **k): return None
    def TextColumn(*a, **k): return None
    def BarColumn(*a, **k): return None
    def TimeElapsedColumn(*a, **k): return None

# --------------------------------------------------------------------------- #
# Cluster defaults (overridable via env/flags) — from the install manifest.
# --------------------------------------------------------------------------- #
CLUSTER = {
    "uma_weights_dir": "/mnt/beegfs/software/mlip-models",
    "uma_default_model": "uma-s-1p1.pt",
    "orca_x86": "/mnt/beegfs/software/orca-6.1.1/x86_64-avx2",
    "orca_arm64": "/mnt/beegfs/software/orca-6.1.1/arm64",
    "openmpi_x86": "/mnt/beegfs/software/openmpi-4.1.8-x86_64",
    "openmpi_arm64": "/mnt/beegfs/software/openmpi-4.1.8-arm64",
    "veloxchem_py": "/mnt/beegfs/software/veloxchem-1.0rc4/x86_64/env/bin/python",
    "gxtb_bin": "/mnt/beegfs/software/g-xtb-2.0.0/x86_64/bin/xtb",
    "metallogen_bin_tmpl": "/mnt/beegfs/software/metallogen-0.0.1/{arch}/env/bin/metallogen",
    "mlip_env_x86": "/mnt/beegfs/software/envs/explorer-x86_64/bin/python",
    "mlip_env_arm64": "/mnt/beegfs/software/envs/explorer-mlip-arm64/bin/python",
    "orca_download": "https://orcaforum.kofo.mpg.de/  (register + download ORCA 6.1.1, license-gated)",
    # GPU envs layer the pyscf GPU venv's site-packages (in-process GPU4PySCF).
    # These mirror the production env-gpu-{x86,arm64} _autode_base_venv.pth recipe.
    "pyscf_venv_x86": "/mnt/beegfs/software/pyscf/venv-x86_64-cuda12",
    "pyscf_venv_arm64": "/mnt/beegfs/software/pyscf/venv-arm64-cuda13",
    "gpu4pyscf_source": "/mnt/beegfs/software/pyscf/gpu4pyscf-source",  # arm64: on the .pth too
    # offline wheelhouse (for no-egress nodes, e.g. GH200) — {arch}-cp{ver} dirs.
    # Staged on BeeGFS OUTSIDE the git repo; override with $AUTODE_WHEELHOUSE.
    "wheelhouse_dir": os.environ.get("AUTODE_WHEELHOUSE", "/mnt/beegfs/software/autode-wheelhouse"),
}
PY_FOR_ARCH = {"x86_64": "3.11", "aarch64": "3.12"}
CORE_DEPS = ["rdkit", "numpy>=1.26,<3", "networkx", "matplotlib", "pillow>=9.5.0",
             "cython", "scipy", "loky", "mendeleev", "nvidia-ml-py3"]


# --------------------------------------------------------------------------- #
# Hardware detection
# --------------------------------------------------------------------------- #
@dataclass
class Hardware:
    arch: str
    glibc: str
    has_gpu: bool
    gpu_names: list = field(default_factory=list)
    cuda: Optional[str] = None
    driver: Optional[str] = None
    online: bool = True   # can we reach PyPI? (GH200 compute nodes cannot)

    @property
    def is_arm(self) -> bool:
        return self.arch == "aarch64"


def _has_egress(host: str = "pypi.org", port: int = 443, timeout: float = 3.0) -> bool:
    import socket
    try:
        socket.create_connection((host, port), timeout=timeout).close()
        return True
    except OSError:
        return False


def detect_hardware() -> Hardware:
    arch = platform.machine()
    try:
        glibc = platform.libc_ver()[1] or "?"
    except Exception:
        glibc = "?"
    online = _has_egress()
    has_gpu, names, cuda, driver = False, [], None, None
    if shutil.which("nvidia-smi"):
        try:
            out = subprocess.run(
                ["nvidia-smi", "--query-gpu=name,driver_version",
                 "--format=csv,noheader"],
                capture_output=True, text=True, timeout=20,
            )
            if out.returncode == 0 and out.stdout.strip():
                has_gpu = True
                for line in out.stdout.strip().splitlines():
                    nm, _, drv = line.partition(",")
                    names.append(nm.strip())
                    driver = drv.strip()
                v = subprocess.run(["nvidia-smi"], capture_output=True, text=True, timeout=20)
                for tok in v.stdout.split():
                    if tok.startswith("12.") or tok.startswith("13.") or tok.startswith("11."):
                        cuda = tok
                        break
        except Exception:
            pass
    return Hardware(arch=arch, glibc=glibc, has_gpu=has_gpu, gpu_names=names,
                    cuda=cuda, driver=driver, online=online)


# --------------------------------------------------------------------------- #
# Checkpoint state (resume/fix)
# --------------------------------------------------------------------------- #
class State:
    def __init__(self, path: Path):
        self.path = path
        self.data = {"completed": {}, "hardware": None, "config": {}}
        if path.exists():
            try:
                self.data = json.loads(path.read_text())
            except Exception:
                console.print(f"[yellow]! ignoring corrupt state file {path}[/]")

    def done(self, step: str) -> bool:
        return step in self.data["completed"]

    def mark(self, step: str, info: dict | None = None):
        self.data["completed"][step] = {"ts": time.strftime("%Y-%m-%dT%H:%M:%S"),
                                        **(info or {})}
        self.save()

    def clear(self, step: str):
        self.data["completed"].pop(step, None)
        self.save()

    def save(self):
        self.path.write_text(json.dumps(self.data, indent=2))


# --------------------------------------------------------------------------- #
# Command runner with live output capture
# --------------------------------------------------------------------------- #
def run(cmd: list, *, env: dict | None = None, cwd: str | None = None,
        timeout: int = 3600, log: Path | None = None) -> tuple[int, str]:
    """Run a command, stream a tail to a log file, return (rc, tail)."""
    full_env = {**os.environ, **(env or {})}
    proc = subprocess.run(cmd, capture_output=True, text=True, env=full_env,
                          cwd=cwd, timeout=timeout)
    out = (proc.stdout or "") + (proc.stderr or "")
    if log:
        log.write_text(out)
    tail = "\n".join(out.strip().splitlines()[-20:])
    return proc.returncode, tail


# --------------------------------------------------------------------------- #
# Installer
# --------------------------------------------------------------------------- #
@dataclass
class Ctx:
    repo: Path
    micromamba: str
    hw: Hardware
    state: State
    prefix: Path
    yes: bool
    editable: bool
    skip_backends: bool
    hf_token: Optional[str]
    logs: Path


class Step:
    def __init__(self, name: str, desc: str, fn: Callable[[Ctx], dict],
                 needed: Callable[[Ctx], bool] = lambda c: True):
        self.name, self.desc, self.fn, self.needed = name, desc, fn, needed


# ---- individual steps ------------------------------------------------------ #
def step_toolchain(c: Ctx) -> dict:
    missing = [t for t in ("cc", "gcc") if not shutil.which(t)]
    if len(missing) == 2:
        raise RuntimeError("No C compiler (cc/gcc) — needed for the Cython "
                           "extensions. Install build-essential / gcc+make.")
    info = {"uv": shutil.which("uv"), "micromamba": c.micromamba,
            "gcc": shutil.which("gcc") or shutil.which("cc"), "make": shutil.which("make")}
    if not info["make"]:
        console.print("[yellow]  ! 'make' not found — most builds are fine without it, "
                      "but install it if a backend fails.[/]")
    return info


def _env_path(c: Ctx) -> Path:
    tag = "gpu-" if c.hw.has_gpu else ""
    arch = "arm64" if c.hw.is_arm else "x86"
    return c.prefix / f"env-{tag}{arch}"


def _pyscf_layer_paths(c: Ctx, pyver: str) -> list[str]:
    """site-packages (+ gpu4pyscf-source on ARM64) to layer into a GPU env,
    mirroring the production env-gpu-* _autode_base_venv.pth recipe."""
    venv = CLUSTER["pyscf_venv_arm64"] if c.hw.is_arm else CLUSTER["pyscf_venv_x86"]
    lines = [f"{venv}/lib64/python{pyver}/site-packages",
             f"{venv}/lib/python{pyver}/site-packages"]
    if c.hw.is_arm:
        lines.append(CLUSTER["gpu4pyscf_source"])   # gpu4pyscf built from source
    return lines


def _wheelhouse(c: Ctx, pyver: str) -> Optional[Path]:
    """Offline wheelhouse dir for this arch/python, if staged."""
    wh = Path(CLUSTER["wheelhouse_dir"]) / f"{c.hw.arch}-cp{pyver.replace('.', '')}"
    return wh if wh.is_dir() and any(wh.iterdir()) else None


def _has_headers(pybin: str) -> bool:
    """True if this interpreter can find its own Python.h (dev headers present).
    The Cython extension won't compile without them."""
    rc, out = run([pybin, "-c",
                   "import sysconfig,os;print(os.path.join(sysconfig.get_path('include'),'Python.h'))"],
                  timeout=30)
    return rc == 0 and Path(out.strip()).exists()


def _pick_base_python(c: Ctx, py: str) -> str:
    """A python{py} interpreter that HAS dev headers (Python.h), to base the
    venv on — else the cconf_gen C-extension can't build. GH200 nodes have
    inconsistent python3.12-dev, so we fall back to a conda/mlip env python
    that reliably bundles headers."""
    cands = [f"/usr/bin/python{py}", shutil.which(f"python{py}")]
    mlip = CLUSTER["mlip_env_arm64"] if c.hw.is_arm else CLUSTER["mlip_env_x86"]
    if Path(mlip).exists():
        rc, out = run([mlip, "-c", "import sys;print('%d.%d'%sys.version_info[:2])"], timeout=30)
        if rc == 0 and out.strip() == py:
            cands.append(mlip)   # conda env python — bundles headers
    for p in cands:
        if p and Path(p).exists() and _has_headers(p):
            return p
    raise RuntimeError(
        f"no python{py} with dev headers (Python.h) found — needed to compile "
        f"the cconf_gen C extension. Install the dev package (e.g. "
        f"`apt install python{py}-dev` / `dnf install python{py}-devel`), or "
        f"ensure a conda env python{py} with headers is available.")


def _mk_venv(c: Ctx, env: Path, py: str) -> None:
    # Use `uv venv` (not `python -m venv`): it needs no system python3-venv /
    # ensurepip package (missing on the GH200 nodes) and works fully offline
    # from the staged uv binary. Base on an interpreter that has dev headers
    # (matching the target ABI) so the Cython extension can compile.
    base = _pick_base_python(c, py)
    rc, tail = run(["uv", "venv", "--python", base, str(env)],
                   log=c.logs / "create_env.log")
    if rc != 0:
        raise RuntimeError(f"uv venv failed:\n{tail}")


def step_create_env(c: Ctx) -> dict:
    env = _env_path(c)
    py = PY_FOR_ARCH.get(c.hw.arch, "3.11")
    epy = env / "bin" / "python"
    if epy.exists() and run([str(epy), "-c", "import sys"], timeout=30)[0] == 0:
        return {"env": str(env), "reused": True, "python": py}
    if env.exists():   # partial/broken env from a failed run — rebuild clean
        shutil.rmtree(env, ignore_errors=True)

    # No egress → we cannot use micromamba (conda-forge download) nor online pip.
    # Fall back to a stdlib venv + the staged offline wheelhouse.
    if not c.hw.online and not c.hw.has_gpu:
        wh = _wheelhouse(c, py)
        if wh is None:
            raise RuntimeError(
                f"this node has no internet egress and no offline wheelhouse for "
                f"{c.hw.arch}-cp{py.replace('.', '')} at {CLUSTER['wheelhouse_dir']}.\n"
                f"Either run on a node with egress, or stage a wheelhouse first "
                f"(pip download ... -d {CLUSTER['wheelhouse_dir']}/{c.hw.arch}-cp{py.replace('.', '')}).")
        _mk_venv(c, env, py)
        return {"env": str(env), "python": py, "kind": "venv-offline", "wheelhouse": str(wh)}

    if c.hw.has_gpu:
        # GPU env: a stdlib venv layered on the pyscf GPU venv's site-packages
        # via a .pth (so in-process GPU4PySCF/torch-cuda import). This is the
        # production recipe — NOT a micromamba env — because the pyscf CUDA
        # stack lives in that venv and --system-site-packages won't chain a venv.
        _mk_venv(c, env, py)
        layer = _pyscf_layer_paths(c, py)
        missing = [p for p in layer if not Path(p).exists()]
        if missing:
            raise RuntimeError(
                "pyscf GPU venv not found — GPU4PySCF layering needs:\n  "
                + "\n  ".join(missing) +
                "\n(install GPU4PySCF first, or run on a node where it is staged.)")
        sp = env / "lib" / f"python{py}" / "site-packages"
        (sp / "_autode_base_venv.pth").write_text("\n".join(layer) + "\n")
        return {"env": str(env), "python": py, "gpu_layer": layer, "kind": "venv"}

    # CPU env: micromamba (conda-forge python + pip).
    rc, tail = run(
        [c.micromamba, "create", "-y", "-p", str(env), "-c", "conda-forge",
         f"python={py}", "pip"],
        env={"MAMBA_ROOT_PREFIX": os.environ.get(
            "MAMBA_ROOT_PREFIX", str(Path(c.micromamba).parent.parent))},
        log=c.logs / "create_env.log",
    )
    if rc != 0:
        raise RuntimeError(f"micromamba create failed:\n{tail}")
    return {"env": str(env), "python": py, "kind": "micromamba"}


def step_install_autode(c: Ctx) -> dict:
    env = _env_path(c)
    py = env / "bin" / "python"
    pyver = PY_FOR_ARCH.get(c.hw.arch, "3.11")
    mode = ["-e"] if c.editable else []

    # Offline (no egress): install deps + build backend from the staged
    # wheelhouse, then build the editable install with --no-build-isolation
    # (build isolation would try to fetch setuptools/cython from PyPI).
    wh = None if c.hw.online else _wheelhouse(c, pyver)
    if wh is not None:
        # uv pip (the env has no pip; uv manages it). --no-index/--find-links
        # → wheelhouse only; --no-build-isolation since setuptools+cython are
        # installed first (build isolation would try to fetch them from PyPI).
        off = ["uv", "pip", "install", "--python", str(py), "--no-index", "--find-links", str(wh)]
        rc, tail = run([*off, "setuptools", "wheel", "cython", *CORE_DEPS],
                       log=c.logs / "install_deps_offline.log")
        if rc != 0:
            raise RuntimeError(f"offline dep install from wheelhouse failed:\n{tail}")
        rc, tail = run([*off, "--no-build-isolation", "--no-deps", *mode, str(c.repo)],
                       log=c.logs / "install_autode.log")
        if rc != 0:
            raise RuntimeError(f"offline autodE editable build failed:\n{tail}")
        return {"editable": bool(mode), "offline": True, "wheelhouse": str(wh)}

    # Online: uv pip is fast; -e recompiles the Cython ext for this py+glibc+arch.
    rc, tail = run(
        ["uv", "pip", "install", "--python", str(py), *mode, str(c.repo)],
        log=c.logs / "install_autode.log",
    )
    if rc != 0:
        # fall back to the env's own pip (no uv) in case of a resolver edge
        rc2, tail2 = run([str(py), "-m", "pip", "install", *mode, str(c.repo)],
                         log=c.logs / "install_autode_pipfallback.log")
        if rc2 != 0:
            raise RuntimeError(f"autodE install failed (uv + pip):\n{tail2}")
    return {"editable": bool(mode), "offline": False}


def step_verify_core(c: Ctx) -> dict:
    py = _env_path(c) / "bin" / "python"
    # cconf_gen is a TOP-LEVEL compiled module (built to the repo root and
    # imported as `from cconf_gen import v` in autode/conformers/conf_gen.py) —
    # NOT autode.conformers.cconf_gen. Verify the real name + the runtime path.
    rc, tail = run([str(py), "-c",
                    "import autode, mendeleev, rdkit, scipy, numpy; "
                    "import cconf_gen; from cconf_gen import v, dvdr; "
                    "from autode.conformers import conf_gen; "
                    "print('autode', autode.__version__)"])
    if rc != 0:
        raise RuntimeError(
            "autodE import / Cython-extension verify failed. The compiled "
            "'cconf_gen' extension may not have built or is not on the path.\n"
            f"{tail}")
    return {"import": tail.strip()}


def step_check_orca(c: Ctx) -> dict:
    orca_dir = os.environ.get("ORCA_DIR") or (
        CLUSTER["orca_arm64"] if c.hw.is_arm else CLUSTER["orca_x86"])
    orca = Path(orca_dir) / "orca"
    ok = orca.exists()
    ver = None
    if ok:
        rc, tail = run([str(orca), "--version"], timeout=30)
        ver = tail.strip().splitlines()[0] if tail.strip() else "present"
    status = "found" if ok else "MISSING"
    color = "green" if ok else "yellow"
    console.print(Panel(
        f"[b]ORCA[/] (external QC backend, license-gated): [{color}]{status}[/]\n"
        f"  ORCA_DIR = {orca_dir}\n"
        + (f"  version  = {ver}\n" if ok else
           f"  → not found. Install ORCA 6.1.1 manually and set ORCA_DIR.\n"
           f"    Download: {CLUSTER['orca_download']}\n"
           f"    OpenMPI 4.1.8 needed on PATH ({CLUSTER['openmpi_arm64' if c.hw.is_arm else 'openmpi_x86']}).\n"),
        title="external: ORCA", border_style=color))
    return {"orca_dir": orca_dir, "found": ok, "version": ver}


def step_backends(c: Ctx) -> dict:
    """Detect + wire the optional method backends the hardware supports."""
    if c.skip_backends:
        return {"skipped": True}
    arch = "arm64" if c.hw.is_arm else "x86_64"
    wired = {}
    # UMA (GPU MLIP) ---------------------------------------------------------
    if c.hw.has_gpu:
        mlip_py = CLUSTER["mlip_env_arm64"] if c.hw.is_arm else CLUSTER["mlip_env_x86"]
        weights = Path(CLUSTER["uma_weights_dir"]) / CLUSTER["uma_default_model"]
        if Path(mlip_py).exists() and weights.exists():
            wired["UMA"] = {"UMA_PROVIDER_PYTHON": mlip_py,
                            "UMA_PROVIDER": str(c.repo / "providers" / "uma_provider.py"),
                            "UMA_MODEL": str(weights)}
        elif not weights.exists():
            console.print("[yellow]  UMA weights not local; off-cluster you'll need an HF token "
                          "(facebook/UMA). Skipping UMA wiring.[/]")
    # GPU4PySCF (in-process; reachable via the create_env pyscf .pth layering) -
    if c.hw.has_gpu:
        env = _env_path(c)
        rc, _ = run([str(env / "bin" / "python"), "-c", "import gpu4pyscf"], timeout=120)
        if rc == 0:
            wired["GPU4PySCF"] = {"(in-process via pyscf layer)": ""}
        else:
            console.print("[yellow]  GPU4PySCF not importable in the env — the pyscf "
                          "GPU-venv layering may be missing (see create_env).[/]")
    # VeloxChem (CPU DFT; glibc>=2.32) --------------------------------------
    vpy = CLUSTER["veloxchem_py"]
    glibc_ok = c.hw.glibc == "?" or tuple(int(x) for x in c.hw.glibc.split(".")[:2]) >= (2, 32)
    if Path(vpy).exists() and (glibc_ok or c.hw.is_arm) and not c.hw.is_arm:
        wired["VeloxChem"] = {"VELOXCHEM_PYTHON": vpy}
    elif Path(vpy).exists() and c.hw.is_arm:
        console.print("[dim]  VeloxChem env is x86-only conda; skip on ARM64.[/]")
    # g-xTB -----------------------------------------------------------------
    if Path(CLUSTER["gxtb_bin"]).exists():
        wired["g-xTB"] = {"GXTB_PATH": CLUSTER["gxtb_bin"]}
    # MetalloGen (metal_builder seed) ---------------------------------------
    mbin = CLUSTER["metallogen_bin_tmpl"].format(arch=arch)
    if Path(mbin).exists():
        wired["MetalloGen"] = {"METALLOGEN_BIN": mbin}
    # report
    t = Table(title="method backends wired for this hardware", show_edge=False)
    t.add_column("backend", style="bold"); t.add_column("env vars")
    for b, kv in wired.items():
        t.add_row(b, ", ".join(kv.keys()))
    if not wired:
        t.add_row("(none)", "core autodE + ORCA only")
    console.print(t)
    return {"wired": wired}


def step_write_activation(c: Ctx) -> dict:
    env = _env_path(c)
    wired = c.state.data["completed"].get("backends", {}).get("wired", {})
    orca_dir = c.state.data["completed"].get("check_orca", {}).get("orca_dir", "")
    lines = ["#!/usr/bin/env bash",
             "# autodE environment — source me. Generated by install.sh.",
             f'export PATH="{env}/bin:$PATH"']
    if orca_dir:
        lines += [f'export ORCA_DIR="{orca_dir}"',
                  f'export PATH="$ORCA_DIR:$PATH"']
        mpi = CLUSTER["openmpi_arm64" if c.hw.is_arm else "openmpi_x86"]
        lines += [f'export PATH="{mpi}/bin:$PATH"',
                  f'export LD_LIBRARY_PATH="{mpi}/lib:${{LD_LIBRARY_PATH:-}}"']
    for _, kv in wired.items():
        for k, v in kv.items():
            lines.append(f'export {k}="{v}"')
    act = env / "activate-autode.sh"
    act.write_text("\n".join(lines) + "\n")
    act.chmod(0o755)
    return {"activation": str(act)}


def step_smoke(c: Ctx) -> dict:
    py = _env_path(c) / "bin" / "python"
    code = ("import autode as ade; from autode.wrappers.UMA import UMA; "
            "from autode.species.metal_builder import get_metal_builder; "
            "print('smoke ok:', ade.__version__, type(get_metal_builder()).__name__)")
    rc, tail = run([str(py), "-c", code])
    if rc != 0:
        console.print(f"[yellow]  smoke test (backends importable) non-fatal warning:\n{tail}[/]")
        return {"smoke": "partial", "detail": tail.strip()[:200]}
    return {"smoke": tail.strip()}


STEPS = [
    Step("toolchain", "verify/prepare build toolchain (gcc, uv, micromamba)", step_toolchain),
    Step("create_env", "create the environment for this hardware (venv+pyscf on GPU, micromamba on CPU)", step_create_env),
    Step("install_autode", "install autodE (uv pip, compiles Cython ext)", step_install_autode),
    Step("verify_core", "verify autodE + Cython import", step_verify_core),
    Step("check_orca", "check external ORCA (path/version) + guidance", step_check_orca),
    Step("backends", "detect + wire method backends for this hardware", step_backends),
    Step("write_activation", "write the activation script (env vars)", step_write_activation),
    Step("smoke", "final smoke test", step_smoke),
]


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def banner(hw: Hardware, prefix: Path, resume: bool):
    g = ", ".join(hw.gpu_names) if hw.gpu_names else "none"
    body = Text.assemble(
        ("autodE installer\n", "bold cyan"),
        (f"arch      {hw.arch}  (glibc {hw.glibc})\n", ""),
        (f"gpu       {g}\n", ""),
        (f"cuda      {hw.cuda or '—'}   driver {hw.driver or '—'}\n", ""),
        (f"egress    {'online' if hw.online else 'OFFLINE → wheelhouse'}\n",
         "" if hw.online else "yellow"),
        (f"prefix    {prefix}\n", ""),
        (f"mode      {'RESUME' if resume else 'fresh'}\n", "yellow" if resume else "green"),
    )
    console.print(Panel(body, border_style="cyan"))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="autodE single-shot installer")
    ap.add_argument("--repo", required=True)
    ap.add_argument("--micromamba", required=True)
    ap.add_argument("--prefix", default=None, help="env install prefix (default: <repo>)")
    ap.add_argument("--resume", action="store_true", help="skip already-completed steps")
    ap.add_argument("--fix", metavar="STEP", help="re-run one step (and everything after)")
    ap.add_argument("--only", metavar="STEP", help="run just one step")
    ap.add_argument("--no-editable", action="store_true", help="install autodE as a copy, not -e")
    ap.add_argument("--skip-backends", action="store_true")
    ap.add_argument("--yes", "-y", action="store_true", help="non-interactive (assume yes)")
    ap.add_argument("--hf-token", default=os.environ.get("HF_TOKEN"))
    a = ap.parse_args(argv)

    repo = Path(a.repo).resolve()
    prefix = Path(a.prefix).resolve() if a.prefix else repo
    logs = prefix / "installer" / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    state = State(prefix / ".autode_install_state.json")

    hw = detect_hardware()
    state.data["hardware"] = asdict(hw)
    state.save()
    banner(hw, prefix, a.resume or bool(a.fix))

    if not (a.yes or a.resume or a.fix or a.only):
        if not Confirm.ask("Proceed with install?", default=True):
            console.print("aborted."); return 1

    ctx = Ctx(repo=repo, micromamba=a.micromamba, hw=hw, state=state, prefix=prefix,
              yes=a.yes, editable=not a.no_editable, skip_backends=a.skip_backends,
              hf_token=a.hf_token, logs=logs)

    # figure out which steps to run
    steps = STEPS
    if a.only:
        steps = [s for s in STEPS if s.name == a.only] or _bad_step(a.only)
    elif a.fix:
        idx = next((i for i, s in enumerate(STEPS) if s.name == a.fix), None)
        if idx is None:
            _bad_step(a.fix)
        for s in STEPS[idx:]:
            state.clear(s.name)
        steps = STEPS

    failed = None
    with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"),
                  BarColumn(), TimeElapsedColumn(), console=console) as prog:
        task = prog.add_task("installing", total=len(steps))
        for s in steps:
            if (a.resume or a.fix) and state.done(s.name) and not (a.only == s.name):
                prog.console.print(f"[dim]• {s.name:16} skipped (done)[/]")
                prog.advance(task); continue
            prog.update(task, description=f"[cyan]{s.name}[/] — {s.desc}")
            try:
                info = s.fn(ctx)
                state.mark(s.name, info)
                prog.console.print(f"[green]✓ {s.name:16}[/] {s.desc}")
            except Exception as e:  # noqa
                prog.console.print(f"[red]✗ {s.name:16} FAILED[/] — {e}")
                prog.console.print(f"[dim]  logs: {logs}[/]")
                failed = s.name
                break
            prog.advance(task)

    console.print(Rule())
    if failed:
        console.print(Panel(
            f"[red]Install stopped at step [b]{failed}[/].[/]\n"
            f"Nothing was lost — fix the issue and re-run:\n"
            f"  [b]./install.sh --resume[/]        continue from here\n"
            f"  [b]./install.sh --fix {failed}[/]  re-run this step (+ later ones)\n"
            f"Logs: {logs}", border_style="red"))
        return 1
    if not a.only:
        env = _env_path(ctx)
        console.print(Panel(
            f"[green b]autodE installed.[/]\n\n"
            f"Activate:  [b]source {env}/activate-autode.sh[/]\n"
            f"Use:       [b]python -c 'import autode as ade'[/]\n"
            f"State:     {state.path}  (delete to reinstall from scratch)",
            title="done", border_style="green"))
    return 0


def _bad_step(name: str):
    console.print(f"[red]unknown step '{name}'. Steps: "
                  f"{', '.join(s.name for s in STEPS)}[/]")
    sys.exit(2)


if __name__ == "__main__":
    raise SystemExit(main())
