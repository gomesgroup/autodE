#!/usr/bin/env bash
# =============================================================================
# autodE single-shot installer — bootstrap
# -----------------------------------------------------------------------------
# Detects hardware, bootstraps modern tooling (uv + micromamba) if missing, then
# hands off to the rich-powered Python driver (installer/autode_install.py) which
# does the real work with progress bars, checkpoint/resume, and backend wiring.
#
#   ./install.sh                 # interactive, auto-detect hardware, install core + supported backends
#   ./install.sh --help          # all options (resume, fix, backend selection, non-interactive, ...)
#   ./install.sh --resume        # continue after a failure (skips completed steps)
#
# Only bash + a POSIX userland is assumed here. Everything else is bootstrapped.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DRIVER="${HERE}/installer/autode_install.py"

c() { printf '\033[%sm' "$1"; }   # tiny colour helper for the pre-rich phase
BOLD="$(c 1)"; DIM="$(c 2)"; GRN="$(c 32)"; YLW="$(c 33)"; RED="$(c 31)"; RST="$(c 0)"
say()  { printf '%s\n' "${BOLD}[autodE-install]${RST} $*"; }
warn() { printf '%s\n' "${YLW}[autodE-install] ! $*${RST}"; }
die()  { printf '%s\n' "${RED}[autodE-install] ✗ $*${RST}" >&2; exit 1; }

[ -f "$DRIVER" ] || die "driver not found at $DRIVER (run from the autodE repo root)"

# --- 1. Hard prerequisites we cannot bootstrap ------------------------------
need_cmd() { command -v "$1" >/dev/null 2>&1; }
MISSING=()
for t in bash uname; do need_cmd "$t" || MISSING+=("$t"); done
# a C toolchain is required for the Cython extensions
if ! need_cmd cc && ! need_cmd gcc; then MISSING+=("gcc/cc (C compiler)"); fi
need_cmd git  || warn "git not found — fine if you already have the source (you do)."
need_cmd curl || need_cmd wget || MISSING+=("curl or wget")
if [ "${#MISSING[@]}" -gt 0 ]; then
  warn "Missing prerequisites the installer can't bootstrap for you:"
  for m in "${MISSING[@]}"; do printf '    - %s\n' "$m"; done
  die "install them (e.g. 'dnf install gcc make' / 'apt install build-essential curl'), then re-run."
fi

fetch() { if need_cmd curl; then curl -LsSf "$1"; else wget -qO- "$1"; fi; }

# --- 2. Locate uv (prefer cluster copy → PATH → bootstrap) ------------------
# Compute/GPU nodes on this cluster have no internet egress, so we prefer the
# shared uv staged on BeeGFS (same pattern as micromamba) and only fall back to
# a network bootstrap when egress is available.
UV=""
_arch="$(uname -m)"
# cluster stages ARM tooling under 'arm64', but uname -m says 'aarch64' — try both
_arch_alt="$_arch"; [ "$_arch" = "aarch64" ] && _arch_alt="arm64"
_runs() { [ -x "$1" ] && "$1" --version >/dev/null 2>&1; }   # exists AND executes here
for cand in \
  "/mnt/beegfs/software/uv/${_arch}/bin/uv" \
  "/mnt/beegfs/software/uv/${_arch_alt}/bin/uv" \
  "$(command -v uv 2>/dev/null || true)" \
  "$HOME/.local/bin/uv" "$HOME/.cargo/bin/uv"; do
  [ -n "$cand" ] && _runs "$cand" && { UV="$cand"; break; }
done
if [ -z "$UV" ]; then
  say "installing ${BOLD}uv${RST} (not found locally or on BeeGFS) ..."
  if fetch https://astral.sh/uv/install.sh | sh >/dev/null 2>&1; then
    export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"
    UV="$(command -v uv 2>/dev/null || true)"
  fi
  [ -n "$UV" ] && [ -x "$UV" ] || die "uv not found and could not be downloaded.
    This node likely has no internet egress. Either run the installer on a node
    with egress (the head node), or point to the shared copy:
      /mnt/beegfs/software/uv/${_arch}/bin/uv"
fi
export PATH="$(dirname "$UV"):$PATH"   # so 'uv' resolves for the exec below
say "uv: ${GRN}$("$UV" --version 2>/dev/null || echo present)${RST} ${DIM}($UV)${RST}"

# --- 3. Locate micromamba (prefer the cluster copy; else bootstrap) ---------
MAMBA="${MICROMAMBA_BIN:-}"
if [ -z "$MAMBA" ]; then
  for cand in \
    "/mnt/beegfs/software/micromamba/${_arch}/bin/micromamba" \
    "/mnt/beegfs/software/micromamba/${_arch_alt}/bin/micromamba" \
    "$(command -v micromamba 2>/dev/null || true)" \
    "$HOME/.local/bin/micromamba"; do
    [ -n "$cand" ] && [ -x "$cand" ] && { MAMBA="$cand"; break; }
  done
fi
if [ -z "$MAMBA" ]; then
  say "installing ${BOLD}micromamba${RST} (not found) ..."
  mkdir -p "$HOME/.local/bin"
  ( cd "$HOME/.local/bin" && fetch "https://micro.mamba.pm/api/micromamba/linux-$(uname -m | sed 's/x86_64/64/;s/aarch64/aarch64/')/latest" | tar -xj bin/micromamba --strip-components=1 ) \
    || die "micromamba bootstrap failed"
  MAMBA="$HOME/.local/bin/micromamba"
fi
say "micromamba: ${GRN}${MAMBA}${RST}"

# --- 4. Hand off to the Python driver ---------------------------------------
# Online: run under `uv run --with rich` for the full rich UI (uv provisions a
# ≥3.9 interpreter + rich in a throwaway env). Offline (no egress — e.g. a
# GH200 compute node): uv can't fetch rich, so run the driver with a system
# python ≥3.9; the driver degrades gracefully to plain output (rich-optional).
# --no-project + cleared VIRTUAL_ENV so uv/py ignore any ambient .venv/pyproject.
unset VIRTUAL_ENV CONDA_PREFIX UV_PROJECT_ENVIRONMENT 2>/dev/null || true

_egress() { timeout 3 bash -c 'exec 3<>/dev/tcp/pypi.org/443' 2>/dev/null; }
if _egress; then
  say "handing off to the installer driver (rich UI) ..."
  exec "$UV" run --quiet --no-project --with "rich>=13" --python ">=3.9" python "$DRIVER" \
    --repo "$HERE" --micromamba "$MAMBA" "$@"
fi

# offline: find a system python >=3.9 to run the driver (rich-optional)
warn "no internet egress — running the installer offline (plain UI, wheelhouse installs)."
PYDRV=""
for p in python3.12 python3.11 python3.10 python3.9 python3; do
  cand="$(command -v "$p" 2>/dev/null || true)"
  [ -n "$cand" ] && "$cand" -c 'import sys;sys.exit(0 if sys.version_info[:2]>=(3,9) else 1)' 2>/dev/null \
    && { PYDRV="$cand"; break; }
done
[ -n "$PYDRV" ] || die "offline and no system python>=3.9 found to run the installer."
say "driver python: ${GRN}${PYDRV}${RST} (offline / rich-optional)"
exec "$PYDRV" "$DRIVER" --repo "$HERE" --micromamba "$MAMBA" "$@"
