"""Pluggable transition-metal (TM) complex *seed* builder for autodE.

RDKit ETKDG has no transition-metal parameters and autodE's own SMILES
``Builder`` assigns coordination geometry from VSEPR/coordination-number alone
(no oxidation state, spin, ligand field, or polyhapticity) and its parser
rejects dative bonds. So neither reliably builds metal-carbonyl / eta^2-alkene
pi-complexes. This module routes a TM SMILES to an external *seed generator*
(MetalloGen by default; molSimplify as fallback) that assembles a chemically
valid starting geometry + explicit connectivity, and re-enters it into a
``Species`` via ``molecule.atoms`` + ``make_graph(bond_list=...)``.

The seed builder does **not** relax the geometry -- it only produces a rough 3D
seed. Relaxation is delegated to whatever autodE ``Method`` the caller selects
(``mol.optimise(method=UMA())`` by default; or GXTB()/VeloxChem()/ORCA()).

Layering (mirrors the organic path, which hands a rough Builder geometry to a
method):

    Molecule(smiles="...[Co]...")
      -> _init_smiles routing fork            (autode/species/molecule.py)
         -> init_metal_smiles(self, smiles)   (this module)
              -> get_metal_builder().build(smiles) -> MetalSeed
                   seed atoms (SMILES order) + bond list + charge/mult
              -> molecule.atoms = seed.atoms
              -> make_graph(molecule, bond_list=seed.bonds)  # keep M-L/eta/dative edges

Adapters are subprocess wrappers (VeloxChem-pattern): they build a command, run
it with a timeout, parse the result, and **fail soft** -- returning a
``MetalSeed(ok=False, error=...)`` rather than raising -- so a node without the
metal stack degrades gracefully instead of crashing a reaction profile.

Environment overrides
---------------------
* ``METALLOGEN_BIN``          -- metallogen console script (arch-detected default)
* ``MOLSIMPLIFY_PYTHON``      -- molSimplify env interpreter (x86_64 only)
* ``MOLSIMPLIFY_BUILD_COMPLEX`` -- molSimplify provider script

References
----------
* Lee et al., "MetalloGen", J. Chem. Inf. Model. 2025 (m-SMILES 3D builder).
* Ioannidis, Gani, Kulik, J. Comput. Chem. 2016 (molSimplify).
"""

from __future__ import annotations

import glob
import json
import os
import platform
import subprocess
import tempfile
from abc import ABC, abstractmethod
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union

from autode.atoms import Atom
from autode.config import Config
from autode.log import logger

# --------------------------------------------------------------------------- #
# Transition-metal predicate (d-block: Sc-Zn, Y-Cd, La, Hf-Hg). Lifted from
# explorer organometallic_geometry.py:107-150.
# --------------------------------------------------------------------------- #
_TRANSITION_METALS = (
    set(range(21, 31)) | set(range(39, 49)) | {57} | set(range(72, 81))
)


def _canonical_smiles(smiles: str) -> Optional[str]:
    """RDKit-canonical SMILES, or None if unparseable (sanitize fallback)."""
    try:
        from rdkit import Chem

        m = Chem.MolFromSmiles(smiles)
        if m is None:
            m = Chem.MolFromSmiles(smiles, sanitize=False)
        return Chem.MolToSmiles(m) if m is not None else None
    except Exception:
        return None


def is_transition_metal_complex(smiles: str) -> bool:
    """True iff ``smiles`` contains a d-block metal **and** has >1 atom.

    A bare metal atom (e.g. ``[Pd]``) is not a "complex" -- it needs no builder
    (a trivial one-atom geometry), so it stays on autodE's normal Builder path.
    """
    try:
        from rdkit import Chem
    except Exception:
        return False

    m = Chem.MolFromSmiles(smiles)
    if m is None:
        m = Chem.MolFromSmiles(smiles, sanitize=False)
    if m is None:
        return False

    has_tm = any(a.GetAtomicNum() in _TRANSITION_METALS for a in m.GetAtoms())
    if not has_tm:
        return False
    try:
        n_atoms = Chem.AddHs(m).GetNumAtoms()
    except Exception:
        n_atoms = m.GetNumAtoms()
    return n_atoms > 1


@dataclass
class MetalSeed:
    """Value object returned by ``MetalComplexBuilder.build``.

    ``atoms``/``bonds`` are consistent: ``bonds`` are index pairs into
    ``atoms``. ``atoms`` are in INPUT-SMILES heavy-atom order (H's appended);
    ``bonds`` include the metal-ligand / dative / eta edges.
    """

    atoms: List[Atom] = field(default_factory=list)
    bonds: List[Tuple[int, int]] = field(default_factory=list)
    charge: int = 0
    mult: int = 1
    ok: bool = True
    error: Optional[str] = None


# --------------------------------------------------------------------------- #
# Registries. Keyed by RDKit-canonical SMILES so lookup is robust to SMILES
# string variation. A general SMILES -> m-SMILES / spec converter is a
# documented follow-on; for now these cover the curated explorer set.
# --------------------------------------------------------------------------- #

# canonical SMILES -> MetalloGen m-SMILES (metal | ligands | geometry).
# Sourced from /mnt/beegfs/software/examples/metallogen/README.md and the
# explorer ORGANOMETALLIC_SPECS set.
_METALLOGEN_RAW: Dict[str, str] = {
    # HCo(CO)3
    "[H][Co]([C-]#[O+])([C-]#[O+])[C-]#[O+]":
        "[Co+1]|[H-:1]|[C-:2]#[O+]|[C-:3]#[O+]|[C-:4]#[O+]|4_tetrahedral",
    # HCo(CO)3(eta2-propene)
    "CC=C->[Co]([C-]#[O+])([C-]#[O+])([C-]#[O+])[H]":
        "[Co+1]|[H-:1]|[C-:2]#[O+]|[C-:3]#[O+]|[C-:4]#[O+]|"
        "C[CH:5]=[CH2:5]|5_trigonal_bipyramidal",
    # propyl-Co(CO)3
    "CCC[Co]([C-]#[O+])([C-]#[O+])[C-]#[O+]":
        "[Co+1]|[CH2-:1]CC|[C-:2]#[O+]|[C-:3]#[O+]|[C-:4]#[O+]|4_tetrahedral",
    # propyl-Co(CO)4
    "CCC[Co]([C-]#[O+])([C-]#[O+])([C-]#[O+])[C-]#[O+]":
        "[Co+1]|[CH2-:1]CC|[C-:2]#[O+]|[C-:3]#[O+]|[C-:4]#[O+]|"
        "[C-:5]#[O+]|5_trigonal_bipyramidal",
    # acyl-Co(CO)3
    "CCCC(=O)[Co]([C-]#[O+])([C-]#[O+])[C-]#[O+]":
        "[Co+1]|CCC[C-:1]=O|[C-:2]#[O+]|[C-:3]#[O+]|[C-:4]#[O+]|4_tetrahedral",
    # acyl-Co(CO)3(H)2
    "CCCC(=O)[Co]([C-]#[O+])([C-]#[O+])([C-]#[O+])([H])[H]":
        "[Co+3]|CCC[C-:1]=O|[C-:2]#[O+]|[C-:3]#[O+]|[C-:4]#[O+]|"
        "[H-:5]|[H-:6]|6_octahedral",
    # Pd(Ph)(I)
    "[Pd](I)c1ccccc1":
        "[Pd+2]|[c-:1]1ccccc1|[I-:2]|2_linear",
    # Pd(H)(I)
    "[Pd](I)[H]":
        "[Pd+2]|[H-:1]|[I-:2]|2_linear",
    # Pd(Ph)(I)(eta2-acrylic)
    "C(=CC(=O)O)->[Pd](I)c1ccccc1":
        "[Pd+2]|[c-:1]1ccccc1|[I-:2]|[CH2:3]=[CH:3]C(=O)O|3_trigonal_planar",
    # Pd(I)(sigma-alkyl)
    "[Pd](I)C(C(=O)O)Cc1ccccc1":
        "[Pd+2]|[I-:1]|[CH-:2](C(=O)O)Cc1ccccc1|2_linear",
}

# canonical SMILES -> molSimplify spec (core / ligands / geometry / oxstate /
# spin). Lifted from explorer organometallic_geometry.py:71-97.
_MOLSIMPLIFY_RAW: Dict[str, dict] = {
    "[H][Co]([C-]#[O+])([C-]#[O+])[C-]#[O+]":
        {"core": "Co", "ligands": "carbonyl,carbonyl,carbonyl,hydride",
         "geometry": "thd", "oxstate": "I", "spin": "1"},
    "[Pd](I)c1ccccc1":
        {"core": "Pd", "ligands": "phenyl,iodide", "geometry": "le",
         "coord": 2, "oxstate": "II", "spin": "1"},
    "[Pd](I)[H]":
        {"core": "Pd", "ligands": "hydride,iodide", "geometry": "le",
         "coord": 2, "oxstate": "II", "spin": "1"},
    "CC=C->[Co]([C-]#[O+])([C-]#[O+])([C-]#[O+])[H]":
        {"core": "Co",
         "ligands": "carbonyl,carbonyl,carbonyl,hydride,CC=C",
         "geometry": "tbp", "coord": 5, "oxstate": "I", "spin": "1"},
    "CCC[Co]([C-]#[O+])([C-]#[O+])[C-]#[O+]":
        {"core": "Co", "ligands": "carbonyl,carbonyl,carbonyl,[CH2]CC",
         "geometry": "thd", "oxstate": "I", "spin": "1"},
    "CCC[Co]([C-]#[O+])([C-]#[O+])([C-]#[O+])[C-]#[O+]":
        {"core": "Co",
         "ligands": "carbonyl,carbonyl,carbonyl,carbonyl,[CH2]CC",
         "geometry": "tbp", "coord": 5, "oxstate": "I", "spin": "1"},
    "CCCC(=O)[Co]([C-]#[O+])([C-]#[O+])[C-]#[O+]":
        # acyl written C-first so molSimplify's connecting atom is the acyl C
        # (kappa1-C acyl). O-first silently gives a Co-O isomer.
        {"core": "Co", "ligands": "carbonyl,carbonyl,carbonyl,[C](CCC)=O",
         "geometry": "thd", "oxstate": "I", "spin": "1"},
    "CCCC(=O)[Co]([C-]#[O+])([C-]#[O+])([C-]#[O+])([H])[H]":
        {"core": "Co",
         "ligands": "carbonyl,carbonyl,carbonyl,hydride,hydride,[C](CCC)=O",
         "geometry": "oct", "coord": 6, "oxstate": "III", "spin": "1"},
    "C(=CC(=O)O)->[Pd](I)c1ccccc1":
        {"core": "Pd", "ligands": "phenyl,iodide,C=CC(=O)O",
         "geometry": "tpl", "coord": 3, "oxstate": "II", "spin": "1"},
    "[Pd](I)C(C(=O)O)Cc1ccccc1":
        {"core": "Pd", "ligands": "iodide,[CH](C(=O)O)Cc1ccccc1",
         "geometry": "le", "coord": 2, "oxstate": "II", "spin": "1"},
}

# Canonicalise keys once at import (robust to SMILES-string variation).
METALLOGEN_MSMILES: Dict[str, str] = {
    (_canonical_smiles(k) or k): v for k, v in _METALLOGEN_RAW.items()
}
MOLSIMPLIFY_SPECS: Dict[str, dict] = {
    (_canonical_smiles(k) or k): v for k, v in _MOLSIMPLIFY_RAW.items()
}


def msmiles_for(smiles: str) -> Optional[str]:
    """Curated MetalloGen m-SMILES for ``smiles`` (canonical match), or None."""
    return METALLOGEN_MSMILES.get(_canonical_smiles(smiles) or smiles)


def spec_for(smiles: str) -> Optional[dict]:
    """Curated molSimplify spec for ``smiles`` (canonical match), or None."""
    return MOLSIMPLIFY_SPECS.get(_canonical_smiles(smiles) or smiles)


# --------------------------------------------------------------------------- #
# Geometry / graph utilities (moved from explorer organometallic_geometry.py:
# 318-440). Backend-agnostic: perceive connectivity from 3D, reorder to a
# target SMILES's heavy-atom order via element-labelled graph (sub)isomorphism.
# --------------------------------------------------------------------------- #


def _parse_xyz(xyz: str):
    """(symbols: list[str], positions: np.ndarray (n,3)) from an xyz string."""
    import numpy as np

    syms: List[str] = []
    pos: List[List[float]] = []
    for line in xyz.splitlines()[2:]:
        p = line.split()
        if len(p) >= 4:
            syms.append(p[0])
            pos.append([float(p[1]), float(p[2]), float(p[3])])
    return syms, np.asarray(pos, dtype=float)


def _rcov(sym: str) -> float:
    from rdkit import Chem

    return Chem.GetPeriodicTable().GetRcovalent(sym)


def _perceive_bonds(symbols, positions, scale: float = 1.3):
    """Explicit bond list (i<j) over ALL atoms; edge if d < scale*(rcov+rcov).

    This is the seed's connectivity. Unlike autodE's default distance-based
    ``make_graph`` (which caps valencies and would drop metal-ligand bonds for
    a hypervalent metal centre), every geometric edge is kept -- so the
    metal-ligand, dative and eta edges survive into the graph.
    """
    import numpy as np

    bonds: List[Tuple[int, int]] = []
    n = len(symbols)
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(positions[i] - positions[j]))
            if d < scale * (_rcov(symbols[i]) + _rcov(symbols[j])):
                bonds.append((i, j))
    return bonds


def _perceive_heavy_graph(symbols, positions, scale: float = 1.3):
    """networkx heavy-atom graph; nodes are original indices with ``el`` attr."""
    import numpy as np
    import networkx as nx

    heavy = [i for i, s in enumerate(symbols) if s != "H"]
    G = nx.Graph()
    for i in heavy:
        G.add_node(i, el=symbols[i])
    for a in range(len(heavy)):
        for b in range(a + 1, len(heavy)):
            i, j = heavy[a], heavy[b]
            d = float(np.linalg.norm(positions[i] - positions[j]))
            if d < scale * (_rcov(symbols[i]) + _rcov(symbols[j])):
                G.add_edge(i, j)
    return G


def _target_heavy_graph(species_smiles: str):
    """networkx heavy-atom graph for a target SMILES, in SMILES atom order.

    Returns (graph, heavy_symbols). Nodes 0..H-1 in the order heavy atoms
    appear in the SMILES (== RDKit atom-iteration order), with an ``el`` attr.
    """
    import networkx as nx
    from rdkit import Chem

    m = Chem.MolFromSmiles(species_smiles)
    if m is None:
        m = Chem.MolFromSmiles(species_smiles, sanitize=False)
    if m is None:
        raise ValueError(
            f"cannot parse target species SMILES: {species_smiles!r}"
        )
    heavy_idx = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() != "H"]
    remap = {idx: k for k, idx in enumerate(heavy_idx)}
    heavy_syms = [m.GetAtomWithIdx(idx).GetSymbol() for idx in heavy_idx]
    G = nx.Graph()
    for k, idx in enumerate(heavy_idx):
        G.add_node(k, el=m.GetAtomWithIdx(idx).GetSymbol())
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in remap and j in remap:
            G.add_edge(remap[i], remap[j])
    return G, heavy_syms


def reorder_species_to_target(symbols, positions, target_smiles: str):
    """Reorder a built species' atoms so heavy atoms follow ``target_smiles``.

    The "reorder keystone". Perceives a heavy-atom graph from the 3D geometry,
    builds the heavy-atom graph of ``target_smiles`` (SMILES order), and finds
    an element-preserving mapping target->built. Heavy atoms are emitted in
    target order; H atoms are appended after. Isomorphism is tried first, then
    a subgraph monomorphism (the perceived graph can carry extra
    metal-coordination edges, e.g. an eta^2-alkene seen as two M-C bonds).

    Returns (symbols_reordered: list[str], positions_reordered: np.ndarray).
    Raises ValueError if no element-preserving mapping exists.
    """
    from networkx.algorithms.isomorphism import (
        GraphMatcher,
        categorical_node_match,
    )

    Gbuilt = _perceive_heavy_graph(symbols, positions)
    Gtgt, tgt_syms = _target_heavy_graph(target_smiles)

    built_heavy = [i for i, s in enumerate(symbols) if s != "H"]
    if len(built_heavy) != Gtgt.number_of_nodes():
        raise ValueError(
            f"heavy-atom count mismatch: built={len(built_heavy)} "
            f"target={Gtgt.number_of_nodes()} for {target_smiles!r}"
        )

    nm = categorical_node_match("el", "X")
    gm = GraphMatcher(Gbuilt, Gtgt, node_match=nm)
    if gm.is_isomorphic():
        b2t = dict(gm.mapping)  # built_idx -> target_pos
    else:
        gm = GraphMatcher(Gbuilt, Gtgt, node_match=nm)
        if gm.subgraph_is_monomorphic():
            b2t = dict(gm.mapping)
        else:
            raise ValueError(
                f"no element-preserving heavy-atom mapping to target "
                f"{target_smiles!r} (built {len(built_heavy)} heavy atoms)"
            )
    t2b = {t: b for b, t in b2t.items()}  # target_pos -> built_idx
    heavy_order = [t2b[t] for t in range(len(tgt_syms))]

    h_indices = [i for i, s in enumerate(symbols) if s == "H"]
    new_order = heavy_order + h_indices

    new_syms = [symbols[i] for i in new_order]
    new_pos = positions[new_order] if len(positions) else positions
    got = [s for s in new_syms if s != "H"]
    if got != tgt_syms:
        raise ValueError(
            f"reorder produced heavy sequence {got} != target {tgt_syms}"
        )
    return new_syms, new_pos


def _target_formula_charge_mult(smiles: str):
    """(element Counter dict, formal charge, multiplicity) for a species SMILES.

    Multiplicity is the parity-minimum consistent with the electron count
    (even -> singlet, odd -> doublet). Returns (None, None, None) if unparseable.
    """
    from rdkit import Chem

    m = Chem.MolFromSmiles(smiles)
    if m is None:
        m = Chem.MolFromSmiles(smiles, sanitize=False)
    if m is None:
        return None, None, None
    try:
        mh = Chem.AddHs(m)
    except Exception:
        mh = m
    counter = dict(Counter(a.GetSymbol() for a in mh.GetAtoms()))
    charge = Chem.GetFormalCharge(m)
    n_e = sum(a.GetAtomicNum() for a in mh.GetAtoms()) - charge
    mult = 1 if n_e % 2 == 0 else 2
    return counter, charge, mult


def _formula_from_symbols(symbols) -> dict:
    return dict(Counter(symbols))


def _atoms_from(symbols, positions) -> List[Atom]:
    return [
        Atom(s, float(p[0]), float(p[1]), float(p[2]))
        for s, p in zip(symbols, positions)
    ]


def _finalise_seed(
    symbols, positions, smiles: str
) -> MetalSeed:
    """Common seed post-processing: charge/mult, formula guard, reorder,
    bond perception, ``Atom`` construction. Backend-agnostic."""
    target, charge, mult = _target_formula_charge_mult(smiles)
    if target is None:
        return MetalSeed(ok=False, error=f"unparseable SMILES: {smiles!r}")

    got = _formula_from_symbols(symbols)
    if got != target:
        return MetalSeed(
            ok=False, charge=charge or 0, mult=mult or 1,
            error=(f"built formula {got} != target {target} for {smiles!r} "
                   "(seed mis-assembly?)"),
        )

    # Reorder heavy atoms to the input-SMILES order (species-level contract).
    # On failure keep the built order (still a valid geometry) with a warning.
    try:
        symbols, positions = reorder_species_to_target(
            symbols, positions, smiles
        )
    except ValueError as e:
        logger.warning(
            f"metal seed reorder-to-SMILES-order failed for {smiles!r} "
            f"({e}); keeping builder order"
        )

    bonds = _perceive_bonds(symbols, positions)
    atoms = _atoms_from(symbols, positions)
    return MetalSeed(atoms=atoms, bonds=bonds, charge=charge, mult=mult, ok=True)


# --------------------------------------------------------------------------- #
# Builders
# --------------------------------------------------------------------------- #


class MetalComplexBuilder(ABC):
    """Seed-builder interface. ``build`` fails soft (never raises)."""

    name: str = "base"

    @property
    @abstractmethod
    def is_available(self) -> bool:
        ...

    @abstractmethod
    def build(self, spec_or_smiles: Union[str, dict]) -> MetalSeed:
        ...


def _arch_tag() -> str:
    return "arm64" if platform.machine() == "aarch64" else "x86_64"


class MetalloGenBuilder(MetalComplexBuilder):
    """Default seed builder: MetalloGen with the UMA MLIP backend (``-c mlip``).

    Arch-portable (x86_64 + ARM64 GH200). Dispatched as a subprocess to the
    metallogen env's console script (so fairchem/torch load in that env, not
    ours). Consumes a curated m-SMILES from ``METALLOGEN_MSMILES``.
    """

    name = "metallogen"

    def __init__(self):
        self._bin = os.environ.get(
            "METALLOGEN_BIN",
            f"/mnt/beegfs/software/metallogen-0.0.1/{_arch_tag()}/env/bin/"
            "metallogen",
        )

    @property
    def is_available(self) -> bool:
        return os.path.exists(self._bin)

    def _run(self, msmiles: str, timeout: int = 900) -> Optional[str]:
        """Run metallogen; return the first result_*.xyz contents, or None."""
        with tempfile.TemporaryDirectory(prefix="metallogen_") as tmp:
            wd = os.path.join(tmp, "scratch")
            sd = os.path.join(tmp, "result")
            os.makedirs(wd, exist_ok=True)
            os.makedirs(sd, exist_ok=True)
            cmd = [
                self._bin, "-s", msmiles, "-c", "mlip",
                "-wd", wd, "-sd", sd, "-r", "1", "-nc", "1",
            ]
            logger.info(f"MetalloGen build: {msmiles}")
            try:
                proc = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=timeout,
                    env=os.environ.copy(),
                )
            except (FileNotFoundError, subprocess.TimeoutExpired) as e:
                logger.error(f"MetalloGen invocation failed: {e}")
                return None
            results = glob.glob(os.path.join(sd, "result_*.xyz"))
            if not results:
                logger.error(
                    "MetalloGen produced no result_*.xyz "
                    f"(rc={proc.returncode}): {proc.stderr.strip()[-400:]}"
                )
                return None

            def _idx(path: str) -> int:
                base = os.path.basename(path)
                try:
                    return int(base[len("result_"):-len(".xyz")])
                except ValueError:
                    return 1 << 30

            with open(sorted(results, key=_idx)[0]) as fh:
                return fh.read()

    def build(self, spec_or_smiles: Union[str, dict]) -> MetalSeed:
        if isinstance(spec_or_smiles, dict):
            return MetalSeed(
                ok=False, error="MetalloGenBuilder needs a SMILES, not a spec"
            )
        smiles = spec_or_smiles
        msmiles = msmiles_for(smiles)
        if msmiles is None:
            return MetalSeed(
                ok=False,
                error=f"no curated MetalloGen m-SMILES for {smiles!r}",
            )
        xyz = self._run(msmiles)
        if not xyz:
            return MetalSeed(ok=False, error="MetalloGen build produced no xyz")
        symbols, positions = _parse_xyz(xyz)
        if not symbols:
            return MetalSeed(ok=False, error="MetalloGen xyz had no atoms")
        return _finalise_seed(symbols, positions, smiles)


class MolSimplifyBuilder(MetalComplexBuilder):
    """Fallback seed builder: molSimplify (x86_64 only). VeloxChem-pattern.

    Lift-and-shift of explorer ``build_metal_complex``: dispatch a JSON spec to
    the molSimplify env's Python, parse the last stdout JSON line. Consumes a
    curated spec from ``MOLSIMPLIFY_SPECS``.
    """

    name = "molsimplify"

    def __init__(self):
        self._py = os.environ.get(
            "MOLSIMPLIFY_PYTHON",
            "/mnt/beegfs/software/molsimplify-2.0.0/x86_64/env/bin/python",
        )
        self._scr = os.environ.get(
            "MOLSIMPLIFY_BUILD_COMPLEX",
            "/mnt/beegfs/software/molsimplify-2.0.0/build_complex.py",
        )

    @property
    def is_available(self) -> bool:
        return os.path.exists(self._py) and os.path.exists(self._scr)

    def _run(self, spec: dict, timeout: int = 600) -> Optional[dict]:
        cmd = [self._py, self._scr, "--json", json.dumps(spec)]
        logger.info(
            f"molSimplify build: core={spec.get('core')} "
            f"geometry={spec.get('geometry')}"
        )
        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=timeout
            )
        except (FileNotFoundError, subprocess.TimeoutExpired) as e:
            logger.error(f"molSimplify invocation failed: {e}")
            return None
        for line in reversed(proc.stdout.strip().splitlines()):
            line = line.strip()
            if line.startswith("{"):
                try:
                    return json.loads(line)
                except json.JSONDecodeError:
                    continue
        logger.error(
            f"molSimplify produced no JSON (rc={proc.returncode}): "
            f"{proc.stderr.strip()[-400:]}"
        )
        return None

    def build(self, spec_or_smiles: Union[str, dict]) -> MetalSeed:
        if isinstance(spec_or_smiles, dict):
            spec, smiles = spec_or_smiles, None
        else:
            smiles = spec_or_smiles
            spec = spec_for(smiles)
            if spec is None:
                return MetalSeed(
                    ok=False,
                    error=f"no curated molSimplify spec for {smiles!r}",
                )
        d = self._run(spec)
        if not d or not d.get("ok") or not d.get("xyz"):
            return MetalSeed(
                ok=False,
                error=(d or {}).get("error", "molSimplify build failed"),
            )
        symbols, positions = _parse_xyz(d["xyz"])
        if not symbols:
            return MetalSeed(ok=False, error="molSimplify xyz had no atoms")
        if smiles is None:
            # Spec-only input: no SMILES to guard/reorder against -- return the
            # raw seed with perceived bonds and neutral charge/mult.
            bonds = _perceive_bonds(symbols, positions)
            return MetalSeed(atoms=_atoms_from(symbols, positions), bonds=bonds)
        return _finalise_seed(symbols, positions, smiles)


# --------------------------------------------------------------------------- #
# Dispatch + re-entry
# --------------------------------------------------------------------------- #

# Instantiated lazily so importing this module never touches the filesystem
# beyond what the constructors read from env.
_ADAPTER_TYPES = {
    "metallogen": MetalloGenBuilder,
    "molsimplify": MolSimplifyBuilder,
}


def get_metal_builder(name: Optional[str] = None) -> MetalComplexBuilder:
    """Resolve a seed builder by name / Config, honouring availability.

    ``name`` (or ``Config.metal_seed_builder``): ``"metallogen"`` (default) |
    ``"molsimplify"`` | ``"auto"`` (arch-portable first, molSimplify fallback).
    Raises ``RuntimeError`` if the requested / any builder is unavailable.
    """
    name = name or getattr(Config, "metal_seed_builder", "metallogen")

    if name == "auto":
        for pref in ("metallogen", "molsimplify"):
            b = _ADAPTER_TYPES[pref]()
            if b.is_available:
                return b
        raise RuntimeError("no metal-complex seed builder available")

    if name not in _ADAPTER_TYPES:
        raise RuntimeError(
            f"unknown metal_seed_builder {name!r}; "
            f"choose from {sorted(_ADAPTER_TYPES)} or 'auto'/'off'"
        )
    return _ADAPTER_TYPES[name]()


def init_metal_smiles(molecule, smiles: str) -> None:
    """Build a TM-complex seed geometry and re-enter it into ``molecule``.

    Sets ``molecule.charge``/``molecule.mult``/``molecule.atoms`` from the seed
    and builds the graph from the seed's explicit bond list (so the
    metal-ligand / dative / eta edges are kept). The geometry is a *seed* --
    the caller relaxes it (``molecule.optimise(method=UMA())``).

    Fails soft: if no builder is available or the build fails, falls back to
    autodE's ``Builder`` path (``init_smiles``) so a node without the metal
    stack still produces a (poorer) geometry rather than crashing.
    """
    from autode.mol_graphs import make_graph
    from autode.smiles.smiles import init_smiles

    try:
        builder = get_metal_builder()
    except RuntimeError as e:
        logger.warning(f"{e}; falling back to autodE Builder for {smiles!r}")
        return init_smiles(molecule, smiles)

    seed = builder.build(smiles)
    if not seed.ok:
        logger.warning(
            f"metal seed build failed for {smiles!r} ({seed.error}); "
            "falling back to autodE Builder"
        )
        return init_smiles(molecule, smiles)

    molecule.charge = seed.charge
    molecule.mult = seed.mult
    molecule.atoms = seed.atoms
    make_graph(molecule, bond_list=seed.bonds)
    logger.info(
        f"Built metal-complex seed for {smiles!r} via {builder.name}: "
        f"{len(seed.atoms)} atoms, {len(seed.bonds)} bonds, "
        f"charge={seed.charge}, mult={seed.mult}"
    )
    return None
