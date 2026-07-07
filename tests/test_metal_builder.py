"""Unit tests for the transition-metal complex seed builder.

The subprocess to MetalloGen / molSimplify is mocked, so these run anywhere
(no GPU, no metal stack). One real end-to-end MetalloGen build is provided
behind a skip guard (set AUTODE_METAL_E2E=1 on a GPU node to run it).
"""

import os
import numpy as np
import pytest

import autode as ade
from autode.config import Config
from autode.species.molecule import Molecule
from autode.species import metal_builder as mb
from autode.species.metal_builder import (
    is_transition_metal_complex,
    get_metal_builder,
    reorder_species_to_target,
    msmiles_for,
    spec_for,
    _target_formula_charge_mult,
    _parse_xyz,
    _finalise_seed,
    MetalSeed,
    MetalloGenBuilder,
    MolSimplifyBuilder,
    METALLOGEN_MSMILES,
    MOLSIMPLIFY_SPECS,
)


# HCo(CO)3 -- canonical registry key and a real (ORCA-XTB2) geometry, in
# MetalloGen's native build order: Co, H, C, O, C, O, C, O.
HCO_SMILES = "[H][Co]([C-]#[O+])([C-]#[O+])[C-]#[O+]"
HCO_XYZ = """8
0\t1\t-22.71085405715
Co -0.41285 -0.066092 -0.062227
H -1.207923 0.799528 0.77992
C 0.517504 -1.130579 -1.110778
O 1.098176 -1.820132 -1.80757
C -0.449258 1.2648 -1.197393
O -0.475038 2.151928 -1.908961
C -0.476846 -1.182991 1.283543
O -0.52852 -1.887565 2.174343
"""


# --------------------------------------------------------------------------- #
# Predicate + registry
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (HCO_SMILES, True),
        ("[Pd](I)c1ccccc1", True),
        ("CC=C->[Co]([C-]#[O+])([C-]#[O+])([C-]#[O+])[H]", True),
        ("CCO", False),
        ("O", False),
        ("[Pd]", False),  # bare metal atom -> not a "complex"
        ("[Na+]", False),  # main-group, single atom
        ("not-a-smiles", False),
    ],
)
def test_is_transition_metal_complex(smiles, expected):
    assert is_transition_metal_complex(smiles) is expected


def test_registry_lookup_canonical():
    assert len(METALLOGEN_MSMILES) == 10
    assert len(MOLSIMPLIFY_SPECS) == 10

    # Lookup works via a non-canonical spelling of the same species.
    m = msmiles_for(HCO_SMILES)
    assert m is not None and m.startswith("[Co+1]") and "4_tetrahedral" in m
    # Same species, different (non-canonical) SMILES ordering resolves too.
    assert msmiles_for("[Co]([H])([C-]#[O+])([C-]#[O+])[C-]#[O+]") == m

    spec = spec_for(HCO_SMILES)
    assert spec is not None and spec["core"] == "Co" and spec["geometry"] == "thd"

    assert msmiles_for("CCO") is None
    assert spec_for("CCO") is None


def test_target_formula_charge_mult():
    formula, charge, mult = _target_formula_charge_mult(HCO_SMILES)
    assert formula == {"Co": 1, "C": 3, "O": 3, "H": 1}
    assert charge == 0
    assert mult == 1


# --------------------------------------------------------------------------- #
# Reorder-to-SMILES-order keystone
# --------------------------------------------------------------------------- #


def test_reorder_to_smiles_order():
    """A scrambled build order comes back with heavy atoms in SMILES order."""
    symbols, positions = _parse_xyz(HCO_XYZ)

    # Scramble: move the metal to the end and shuffle a CO pair.
    order = [2, 3, 4, 5, 1, 6, 7, 0]  # C O C O H C O Co
    sc_syms = [symbols[i] for i in order]
    sc_pos = positions[order]

    new_syms, new_pos = reorder_species_to_target(sc_syms, sc_pos, HCO_SMILES)

    heavy = [s for s in new_syms if s != "H"]
    # SMILES heavy order for [H][Co](CO)(CO)CO is Co, C, O, C, O, C, O.
    assert heavy == ["Co", "C", "O", "C", "O", "C", "O"]
    # H is appended after the heavy atoms.
    assert new_syms[-1] == "H"
    assert new_pos.shape == (8, 3)


def test_reorder_raises_on_mismatch():
    symbols, positions = _parse_xyz(HCO_XYZ)
    with pytest.raises(ValueError):
        reorder_species_to_target(symbols, positions, "CCO")  # wrong target


# --------------------------------------------------------------------------- #
# _finalise_seed (formula guard, reorder, bond perception)
# --------------------------------------------------------------------------- #


def test_finalise_seed_bonds_and_order():
    symbols, positions = _parse_xyz(HCO_XYZ)
    seed = _finalise_seed(symbols, positions, HCO_SMILES)

    assert seed.ok
    assert seed.charge == 0 and seed.mult == 1
    assert len(seed.atoms) == 8
    # Heavy atoms first in SMILES order; H last.
    assert [a.label for a in seed.atoms] == [
        "Co", "C", "O", "C", "O", "C", "O", "H"
    ]

    # 3 Co-C + 3 C-O + 1 Co-H = 7 edges.
    assert len(seed.bonds) == 7
    edges = {frozenset(b) for b in seed.bonds}
    # Co is atom 0; it must bond the three carbonyl carbons (1,3,5) and the
    # hydride (7) -- the M-L / dative edges survive.
    assert {frozenset((0, 1)), frozenset((0, 3)),
            frozenset((0, 5)), frozenset((0, 7))} <= edges
    # Each carbonyl C-O.
    assert {frozenset((1, 2)), frozenset((3, 4)),
            frozenset((5, 6))} <= edges


def test_finalise_seed_formula_guard():
    """A mis-assembled geometry (wrong formula) is rejected."""
    symbols, positions = _parse_xyz(HCO_XYZ)
    # Drop the hydride -> C3O3Co != target C3O3CoH.
    symbols2 = symbols[:1] + symbols[2:]
    positions2 = np.vstack([positions[:1], positions[2:]])
    seed = _finalise_seed(symbols2, positions2, HCO_SMILES)
    assert not seed.ok
    assert "formula" in (seed.error or "")


# --------------------------------------------------------------------------- #
# MetalloGen adapter (subprocess mocked)
# --------------------------------------------------------------------------- #


def _mock_metallogen_run(sd_writes_xyz=HCO_XYZ):
    """Return a subprocess.run replacement that writes result_1.xyz into -sd."""

    class _Proc:
        returncode = 0
        stdout = ""
        stderr = ""

    def _run(cmd, *a, **k):
        sd = cmd[cmd.index("-sd") + 1]
        os.makedirs(sd, exist_ok=True)
        with open(os.path.join(sd, "result_1.xyz"), "w") as fh:
            fh.write(sd_writes_xyz)
        return _Proc()

    return _run


def test_metallogen_build_mocked(monkeypatch):
    monkeypatch.setattr(mb.subprocess, "run", _mock_metallogen_run())

    b = MetalloGenBuilder()
    monkeypatch.setattr(b, "_bin", "/does/not/matter")  # bypass is_available
    seed = b.build(HCO_SMILES)

    assert seed.ok
    assert len(seed.atoms) == 8
    assert seed.charge == 0 and seed.mult == 1
    assert [a.label for a in seed.atoms][0] == "Co"
    assert len(seed.bonds) == 7


def test_metallogen_build_no_registry_entry():
    seed = MetalloGenBuilder().build("[Fe](Cl)(Cl)Cl")  # no curated m-SMILES
    assert not seed.ok
    assert "m-SMILES" in (seed.error or "")


def test_metallogen_build_rejects_spec_dict():
    seed = MetalloGenBuilder().build({"core": "Co"})
    assert not seed.ok


# --------------------------------------------------------------------------- #
# molSimplify adapter (subprocess mocked)
# --------------------------------------------------------------------------- #


def test_molsimplify_build_mocked(monkeypatch):
    import json

    class _Proc:
        returncode = 0
        stderr = ""
        stdout = "chatter\n" + json.dumps(
            {"ok": True, "xyz": HCO_XYZ, "natoms": 8}
        )

    monkeypatch.setattr(mb.subprocess, "run", lambda *a, **k: _Proc())

    seed = MolSimplifyBuilder().build(HCO_SMILES)
    assert seed.ok
    assert len(seed.atoms) == 8
    assert seed.charge == 0 and seed.mult == 1
    assert len(seed.bonds) == 7


def test_molsimplify_build_soft_fail(monkeypatch):
    class _Proc:
        returncode = 1
        stderr = "boom"
        stdout = ""

    monkeypatch.setattr(mb.subprocess, "run", lambda *a, **k: _Proc())
    seed = MolSimplifyBuilder().build(HCO_SMILES)
    assert not seed.ok


# --------------------------------------------------------------------------- #
# Dispatch
# --------------------------------------------------------------------------- #


def test_get_metal_builder_by_name():
    assert get_metal_builder("metallogen").name == "metallogen"
    assert get_metal_builder("molsimplify").name == "molsimplify"
    with pytest.raises(RuntimeError):
        get_metal_builder("nonsense")


def test_get_metal_builder_default_from_config(monkeypatch):
    monkeypatch.setattr(Config, "metal_seed_builder", "molsimplify")
    assert get_metal_builder().name == "molsimplify"


# --------------------------------------------------------------------------- #
# Re-entry contract + routing (via Molecule)
# --------------------------------------------------------------------------- #


class _FakeBuilder:
    name = "fake"

    def __init__(self, seed):
        self._seed = seed

    @property
    def is_available(self):
        return True

    def build(self, smiles):
        return self._seed


def test_reentry_contract(monkeypatch):
    """A built seed lands on the Molecule: atoms, charge, mult, graph edges."""
    symbols, positions = _parse_xyz(HCO_XYZ)
    seed = _finalise_seed(symbols, positions, HCO_SMILES)
    monkeypatch.setattr(
        mb, "get_metal_builder", lambda name=None: _FakeBuilder(seed)
    )

    mol = Molecule(smiles=HCO_SMILES)

    assert mol.n_atoms == 8
    assert mol.charge == 0 and mol.mult == 1
    assert mol.atoms[0].label == "Co"
    # Graph carries exactly the seed's explicit bonds (incl. dative/eta edges).
    assert mol.graph.number_of_edges() == len(seed.bonds)
    # Co node has degree 4 (3 CO + hydride) -- M-L edges were NOT valence-capped.
    assert mol.graph.degree[0] == 4


def test_routing_calls_metal_path(monkeypatch):
    called = {}
    monkeypatch.setattr(
        "autode.species.molecule.init_metal_smiles",
        lambda self, smiles: called.setdefault("metal", smiles),
    )
    Molecule(smiles=HCO_SMILES)
    assert called.get("metal") == HCO_SMILES


def test_routing_off_uses_builder(monkeypatch):
    monkeypatch.setattr(Config, "metal_seed_builder", "off")
    called = {}
    monkeypatch.setattr(
        "autode.species.molecule.init_metal_smiles",
        lambda self, smiles: called.setdefault("metal", smiles),
    )
    monkeypatch.setattr(
        "autode.species.molecule.init_smiles",
        lambda self, smiles: called.setdefault("builder", smiles),
    )
    Molecule(smiles=HCO_SMILES)
    assert "metal" not in called
    assert called.get("builder") == HCO_SMILES


def test_routing_organic_untouched(monkeypatch):
    called = {}
    monkeypatch.setattr(
        "autode.species.molecule.init_metal_smiles",
        lambda self, smiles: called.setdefault("metal", smiles),
    )
    mol = Molecule(smiles="CCO")  # ethanol -> organic RDKit path
    assert "metal" not in called
    assert mol.n_atoms == 9  # C2H6O


def test_init_metal_smiles_soft_fallback(monkeypatch):
    """If the builder fails, we fall back to autodE's Builder (no crash)."""
    monkeypatch.setattr(
        mb, "get_metal_builder",
        lambda name=None: _FakeBuilder(MetalSeed(ok=False, error="nope")),
    )
    fell_back = {}
    monkeypatch.setattr(
        "autode.smiles.smiles.init_smiles",
        lambda self, smiles: fell_back.setdefault("builder", smiles),
    )
    from autode.species.metal_builder import init_metal_smiles

    m = Molecule(name="x")
    init_metal_smiles(m, HCO_SMILES)
    assert fell_back.get("builder") == HCO_SMILES


# --------------------------------------------------------------------------- #
# Real end-to-end (skipped unless explicitly enabled on a GPU node)
# --------------------------------------------------------------------------- #

_RUN_E2E = os.environ.get("AUTODE_METAL_E2E") == "1"


@pytest.mark.skipif(
    not (_RUN_E2E and MetalloGenBuilder().is_available),
    reason="set AUTODE_METAL_E2E=1 on a node with MetalloGen (-c mlip needs a GPU)",
)
def test_metallogen_real_e2e():
    """Build HCo(CO)3 through the real MetalloGen -c mlip subprocess."""
    seed = MetalloGenBuilder().build(HCO_SMILES)
    assert seed.ok, seed.error
    assert len(seed.atoms) == 8
    assert seed.charge == 0 and seed.mult == 1
    assert seed.atoms[0].label == "Co"

    coords = np.array([a.coord for a in seed.atoms])
    # Co-C distances (carbonyl carbons are atoms 1,3,5 after reorder) ~1.7-1.9 A.
    co_c = [np.linalg.norm(coords[0] - coords[i]) for i in (1, 3, 5)]
    assert all(1.6 < d < 2.0 for d in co_c), co_c
    # Bond list has the 4 M-L edges from the metal.
    assert sum(1 for b in seed.bonds if 0 in b) == 4
