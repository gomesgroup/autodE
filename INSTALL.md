# Installing autodE (single-shot, hardware-aware)

One command sets up autodE correctly for whatever hardware you're on — CPU or GPU,
x86_64 or ARM64 — bootstrapping modern tooling, compiling the Cython extensions,
and wiring the optional method backends it detects.

```bash
git clone https://github.com/gomesgroup/autode.git
cd autode
./install.sh
```

That's it. The installer detects your hardware, creates an isolated env, installs
autodE + its dependencies, checks for external tools (ORCA), wires the MLIP/QC
backends your hardware supports, and writes an activation script. Then:

```bash
source env-<hw>/activate-autode.sh        # path printed at the end
python -c "import autode as ade; print(ade.__version__)"
```

## What it does (and why one script is enough)

`install.sh` is a thin bootstrap; the real work is a `rich`-powered Python driver
(`installer/autode_install.py`) run under `uv run --with rich` — so you get
progress bars and clean logging with zero global installs.

| step | what |
|------|------|
| `toolchain` | verifies a C compiler (for the Cython ext), locates **uv** + **micromamba** (prefers the shared copies on BeeGFS; bootstraps only if missing + egress) |
| `create_env` | **CPU:** `micromamba create` (Python 3.11 x86 / 3.12 ARM64). **GPU:** a stdlib `venv` layered on the pyscf CUDA venv's `site-packages` via a `.pth` — so in-process **GPU4PySCF**/torch-cuda import (the production skala recipe) |
| `install_autode` | `uv pip install -e` autodE — **recompiles the Cython `cconf_gen` extension** for this exact Python+glibc+arch (this is how boltzmann's glibc-2.28 build stays separate from the Ubuntu glibc-2.39 nodes) |
| `verify_core` | imports autodE + the compiled `cconf_gen` to prove the build |
| `check_orca` | checks the external **ORCA** path/version and prints manual-install guidance if missing |
| `backends` | detects + wires the MLIP/QC backends your hardware supports |
| `write_activation` | writes `env-<hw>/activate-autode.sh` (PATH + all backend env vars) |
| `smoke` | final import test incl. the wired backends |

## Prerequisites (the only things it can't bootstrap)

- **bash**, a **C compiler** (`gcc`/`cc`) + ideally `make` — for the Cython extensions.
  - Rocky/RHEL: `sudo dnf install gcc gcc-c++ make`
  - Debian/Ubuntu: `sudo apt install build-essential`
- **git**.

`uv` and `micromamba` are already staged on BeeGFS (`/mnt/beegfs/software/{uv,micromamba}/`)
and are found automatically on every node — no download needed. The env and all Python
deps are installed for you.

## Where to run it (egress + hardware must match the target env)

Package downloads need internet egress, and the Cython extension is compiled for the
**node you run on** — so run the installer **on a node that both has egress and is the
architecture/GPU you're building for**:

| target env | run the installer on | why |
|------------|----------------------|-----|
| x86 CPU (`env-x86`) | **gpg-head** | x86, glibc 2.39, egress; the `.so` works on all Ubuntu x86 nodes |
| x86 GPU (`env-gpu-x86`) | **gpg-boltzmann** (`srun -p gpu-x86 --gres=gpu:1`) | its glibc-2.28 `.so` + A100/A30; layers the x86 pyscf venv |
| ARM64 GPU (`env-gpu-arm64`) | a **GH200** (`srun -p gpu-gh --gres=gpu:1`) | aarch64 `.so`; layers the ARM64 pyscf venv |

Compute nodes with **no egress** (the EPYC `cpu-epyc` nodes) can't download packages — the
installer detects a missing `uv`/egress and tells you to use an egress node instead.

## Hardware scenarios (auto-detected)

| you're on | env | backends wired |
|-----------|-----|----------------|
| x86_64, no GPU (EPYC / head) | `env-x86` | ORCA*, VeloxChem, g-xTB |
| x86_64 + GPU (boltzmann) | `env-gpu-x86` | + UMA (MLIP), GPU4PySCF, MetalloGen `-c mlip` |
| ARM64 + GPU (GH200) | `env-gpu-arm64` | + UMA, GPU4PySCF, MetalloGen, VeloxChem(CPU) |

\* ORCA is external/license-gated — see below.

## Optional method backends

autodE's wrapper methods shell out to external tools; the installer wires the ones
present on this cluster via env vars written into the activation script:

| backend | GPU | wired via | source |
|---------|-----|-----------|--------|
| **ORCA** (DFT/HF/xTB) | no | `ORCA_DIR` + OpenMPI on PATH | **install manually** (below) |
| **UMA** (MLIP) | yes | `UMA_PROVIDER_PYTHON`, `UMA_PROVIDER`, `UMA_MODEL` | fairchem env + local weights `/mnt/beegfs/software/mlip-models/` |
| **GPU4PySCF** | yes | GPU env overlay | pyscf GPU venv |
| **VeloxChem** (CPU DFT) | no | `VELOXCHEM_PYTHON` | needs glibc ≥ 2.32 (not boltzmann) |
| **g-xTB** | no | `GXTB_PATH` | g-xtb-2.0.0 binary |
| **MetalloGen** (metal builder) | yes (`-c mlip`) | `METALLOGEN_BIN` | metallogen env |

Select with `method=UMA()` / `Config.lcode="uma"` etc.; the metal_builder uses
MetalloGen by default (`Config.metal_seed_builder`).

### ORCA (manual, license-gated)

ORCA can't be auto-downloaded (registration + license). Install it, then set `ORCA_DIR`:

1. Register + download **ORCA 6.1.1** from <https://orcaforum.kofo.mpg.de/>.
2. On this cluster it's already at `/mnt/beegfs/software/orca-6.1.1/{x86_64-avx2,arm64}` with
   OpenMPI 4.1.8 (`/mnt/beegfs/software/openmpi-4.1.8-{x86_64,arm64}`) — the installer finds it.
3. Elsewhere: `export ORCA_DIR=/path/to/orca` before running, or edit the activation script.

The installer **checks the ORCA path + version at install time** and prints guidance if missing —
it never blocks the rest of the install.

### UMA weights / HF token

UMA model weights are **local** on the cluster (`/mnt/beegfs/software/mlip-models/uma-s-1p1.pt`),
so no download is needed here. Off-cluster, `facebook/UMA` is HF-gated — pass `--hf-token <token>`
(or `export HF_TOKEN=...`) and the installer will fetch them.

## Resume / fix (never restart from zero)

Every step is checkpointed to `<prefix>/.autode_install_state.json`. If anything fails:

```bash
./install.sh --resume          # continue from the failed step
./install.sh --fix create_env  # re-run one step (and everything after it)
./install.sh --only backends    # run just one step
```

Delete `.autode_install_state.json` to force a clean reinstall.

## Options

```
--prefix DIR      install the env somewhere other than the repo
--no-editable     install autodE as a copy (default is editable -e)
--skip-backends   core autodE + ORCA only
--yes / -y        non-interactive
--hf-token TOK    HF token for off-cluster UMA weights
--resume / --fix STEP / --only STEP
```

## Verifying

```bash
source env-<hw>/activate-autode.sh
python -c "import autode as ade; from autode.species.metal_builder import get_metal_builder; \
          print(ade.__version__, type(get_metal_builder()).__name__)"
```
