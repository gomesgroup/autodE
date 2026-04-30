"""
g-xTB driver for autodE.

``g-xTB`` (Grimme et al.) is a semiempirical method that approximates
ωB97M-V/def2-TZVPPD properties for elements Z = 1–103. It is invoked through
the same ``xtb`` binary as standard GFN1/2-xTB but with an additional
``--gxtb`` flag.

Upstream autodE 1.4.5 has no first-class g-xTB support, so this module
provides a thin subclass of :class:`autode.wrappers.XTB.XTB` that:

1. points at the g-xTB 2.0.0 binary by default
   (``/mnt/beegfs/software/g-xtb-2.0.0/x86_64/bin/xtb`` on the gpg cluster,
   overridable via the ``GXTB_PATH`` environment variable), and
2. always prepends ``--gxtb`` to the runtime flag list, while suppressing
   the default ``--gfn N`` flag (g-xTB has its own parameter set).

Architecture
------------

On x86_64 nodes this is a native binary. On ARM64 (Grace Hopper) nodes the
x86_64 binary runs under QEMU ``binfmt_misc`` because Grimme does not ship
a Linux ARM64 build; energies match x86_64 native to ~1e-11 Eh.

Usage
-----

>>> import autode as ade
>>> from autode.wrappers.GXTB import GXTB
>>>
>>> rxn = ade.Reaction("CC(=O)O.CO>>CC(=O)OC.O", name="esterification")
>>> rxn.calculate_reaction_profile(
...     lmethod=GXTB(),
...     hmethod=ade.methods.ORCA(),
... )

Or for a single optimization:

>>> mol = ade.Molecule(smiles="CCO")
>>> mol.optimise(method=GXTB())

References
----------

* g-xTB preprint: https://chemrxiv.org/engage/chemrxiv/article-details/685434533ba0887c335fc974
* g-xTB releases: https://github.com/grimme-lab/g-xtb
"""

import os

from autode.wrappers.XTB import XTB
from autode.wrappers.keywords import OptKeywords, GradientKeywords

# Default cluster install path; overridable via GXTB_PATH env var so the same
# class works inside Docker images that bind-mount g-xTB at a different path.
_DEFAULT_GXTB_PATH = "/mnt/beegfs/software/g-xtb-2.0.0/x86_64/bin/xtb"


class GXTB(XTB):
    """g-xTB 2.0.0 driver — subclass of :class:`XTB` with ``--gxtb`` injected."""

    def __init__(self):
        super().__init__()

        # Distinct name so Config.lcode = "gxtb" resolves to this class
        # (and not to plain XTB). The on-disk executable is still "xtb".
        self._name = "gxtb"

        # Override the binary path: use g-xTB-bundled xtb (which also handles
        # vanilla --gfn 0/1/2 identically to standard xtb 6.7.1, so this is a
        # safe upgrade). Honour GXTB_PATH for Docker / non-cluster use.
        self.path = os.environ.get("GXTB_PATH") or _DEFAULT_GXTB_PATH

        # g-xTB has its own parameter set baked in; --gfn N would be ignored
        # at best and incompatible at worst. Force it off.
        self.gfn_version = None

        # Cite g-xTB, not standard xTB.
        self.doi_list = [
            # Froitzheim, Grimme et al., "g-xTB: a general-purpose extended
            # tight-binding electronic structure method for the elements
            # H to Lr (Z = 1–103)", ChemRxiv (2025).
            "10.26434/chemrxiv-2025-bjxvt",
        ]

    def __repr__(self):
        return f"GXTB(available = {self.is_available})"

    def execute(self, calc):
        """Prepend ``--gxtb`` to keywords, then delegate to ``XTB.execute()``.

        Subtlety: :meth:`XTB.execute` adds ``--opt`` / ``--grad`` implicitly
        only when ``calc.input.keywords`` is empty. Naively prepending
        ``--gxtb`` would make the list non-empty and silently turn
        optimizations into single-points. We therefore add the implicit
        ``--opt`` / ``--grad`` ourselves when the original list was empty.
        """
        kw = calc.input.keywords
        existing = list(kw)
        if "--gxtb" in existing:
            return super().execute(calc)

        # Reproduce XTB.execute()'s implicit-flag logic that we'd otherwise suppress.
        extra: list[str] = []
        if len(existing) == 0:
            if isinstance(kw, OptKeywords):
                extra.append("--opt")
            elif isinstance(kw, GradientKeywords):
                extra.append("--grad")

        kw_cls = type(kw)
        new_kw = kw_cls(["--gxtb"] + extra + existing)
        # Preserve max_opt_cycles for OptKeywords (None-safe).
        if hasattr(kw, "max_opt_cycles"):
            try:
                new_kw.max_opt_cycles = kw.max_opt_cycles
            except Exception:
                pass
        calc.input.keywords = new_kw
        return super().execute(calc)
