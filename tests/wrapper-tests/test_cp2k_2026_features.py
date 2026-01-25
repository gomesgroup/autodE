#!/usr/bin/env python3
"""
Unit tests for CP2K 2026.1 features in the autodE wrapper.

Tests the new functionality without requiring the full autodE environment
or an actual CP2K installation.
"""
import sys
import os
import tempfile

# Add autode to path
sys.path.insert(0, '/mnt/beegfs/software/autode')

# Mock the autode imports we need
class MockConfig:
    class CP2K:
        path = "/usr/bin/cp2k.psmp"
        keywords = None
        implicit_solvation_type = None

# Patch Config before importing CP2K wrapper
import autode.config
autode.config.Config = MockConfig

# Now we can test the wrapper logic
from autode.wrappers.CP2K import CP2K


def test_aug_molopt_basis_mapping():
    """Test AUG_MOLOPT basis set mapping (CP2K 2026.1+)"""
    print("Testing AUG_MOLOPT basis set mapping...")

    cp2k = CP2K()

    # Test augmented basis sets
    assert cp2k._cp2k_basis("aug-dzvp", "C") == "AUG-DZVP-MOLOPT-GTH"
    assert cp2k._cp2k_basis("aug-tzvp", "O") == "AUG-TZVP-MOLOPT-GTH"
    assert cp2k._cp2k_basis("AUG-DZVP", "H") == "AUG-DZVP-MOLOPT-GTH"
    assert cp2k._cp2k_basis("aug", "N") == "AUG-DZVP-MOLOPT-GTH"  # Default augmented

    # Test standard basis sets still work
    assert cp2k._cp2k_basis("dzvp", "C") == "DZVP-MOLOPT-GTH"
    assert cp2k._cp2k_basis("tzvp", "O") == "TZVP-MOLOPT-GTH"
    assert cp2k._cp2k_basis("tzv2p", "H") == "TZV2P-MOLOPT-GTH"
    assert cp2k._cp2k_basis("unknown", "C") == "DZVP-MOLOPT-GTH"  # Default

    print("  ✓ AUG_MOLOPT basis set mapping works correctly")


def test_sftddft_detection():
    """Test SF-TDDFT keyword detection (CP2K 2026.1+)"""
    print("Testing SF-TDDFT keyword detection...")

    cp2k = CP2K()

    # Mock calculation executor with keywords
    class MockKeyword:
        def __init__(self, kw):
            self._kw = kw
        def __str__(self):
            return self._kw

    class MockInput:
        def __init__(self, keywords):
            self.keywords = [MockKeyword(k) for k in keywords]

    class MockCalc:
        def __init__(self, keywords):
            self.input = MockInput(keywords)

    # Test SF-TDDFT detection
    assert cp2k._has_sftddft(MockCalc(["pbe", "sf-tddft"])) == True
    assert cp2k._has_sftddft(MockCalc(["pbe", "sftddft"])) == True
    assert cp2k._has_sftddft(MockCalc(["pbe", "spinflip"])) == True
    assert cp2k._has_sftddft(MockCalc(["pbe", "tddft"])) == False
    assert cp2k._has_sftddft(MockCalc(["pbe", "opt"])) == False

    print("  ✓ SF-TDDFT keyword detection works correctly")


def test_tddft_detection():
    """Test standard TDDFT keyword detection"""
    print("Testing TDDFT keyword detection...")

    cp2k = CP2K()

    class MockKeyword:
        def __init__(self, kw):
            self._kw = kw
        def __str__(self):
            return self._kw

    class MockInput:
        def __init__(self, keywords):
            self.keywords = [MockKeyword(k) for k in keywords]

    class MockCalc:
        def __init__(self, keywords):
            self.input = MockInput(keywords)

    # Test TDDFT detection
    assert cp2k._has_tddft(MockCalc(["pbe", "tddft"])) == True
    assert cp2k._has_tddft(MockCalc(["pbe", "tddfpt"])) == True
    assert cp2k._has_tddft(MockCalc(["pbe", "opt"])) == False

    print("  ✓ TDDFT keyword detection works correctly")


def test_nstates_parsing():
    """Test excited states count parsing"""
    print("Testing nstates parsing...")

    cp2k = CP2K()

    class MockKeyword:
        def __init__(self, kw):
            self._kw = kw
        def __str__(self):
            return self._kw

    class MockInput:
        def __init__(self, keywords):
            self.keywords = [MockKeyword(k) for k in keywords]

    class MockCalc:
        def __init__(self, keywords):
            self.input = MockInput(keywords)

    # Test nstates parsing
    assert cp2k._get_nstates(MockCalc(["pbe", "tddft", "nstates=10"])) == 10
    assert cp2k._get_nstates(MockCalc(["pbe", "tddft", "nroots=15"])) == 15
    assert cp2k._get_nstates(MockCalc(["pbe", "tddft"])) == 5  # Default

    print("  ✓ nstates parsing works correctly")


def test_cg_optimizer_detection():
    """Test CG optimizer detection (for 3PNT linesearch fix)"""
    print("Testing CG optimizer detection...")

    cp2k = CP2K()

    class MockKeyword:
        def __init__(self, kw):
            self._kw = kw
        def __str__(self):
            return self._kw

    class MockInput:
        def __init__(self, keywords):
            self.keywords = [MockKeyword(k) for k in keywords]

    class MockCalc:
        def __init__(self, keywords):
            self.input = MockInput(keywords)

    # Test CG optimizer detection
    assert cp2k._has_cg_optimizer(MockCalc(["pbe", "cg", "opt"])) == True
    assert cp2k._has_cg_optimizer(MockCalc(["pbe", "conjugate", "opt"])) == True
    assert cp2k._has_cg_optimizer(MockCalc(["pbe", "opt"])) == False
    assert cp2k._has_cg_optimizer(MockCalc(["pbe", "bfgs", "opt"])) == False

    print("  ✓ CG optimizer detection works correctly")


def test_run_type_detection():
    """Test run type detection including TDDFT"""
    print("Testing run type detection...")

    cp2k = CP2K()

    class MockKeyword:
        def __init__(self, kw):
            self._kw = kw
        def __str__(self):
            return self._kw

    class MockInput:
        def __init__(self, keywords):
            self.keywords = [MockKeyword(k) for k in keywords]

    class MockCalc:
        def __init__(self, keywords):
            self.input = MockInput(keywords)

    # Test run type detection
    assert cp2k._get_run_type(MockCalc(["pbe", "opt"])) == "GEO_OPT"
    assert cp2k._get_run_type(MockCalc(["pbe", "ts"])) == "TRANSITION_STATE"
    assert cp2k._get_run_type(MockCalc(["pbe", "force"])) == "ENERGY_FORCE"
    assert cp2k._get_run_type(MockCalc(["pbe", "freq"])) == "VIBRATIONAL_ANALYSIS"
    assert cp2k._get_run_type(MockCalc(["pbe", "sf-tddft"])) == "ENERGY"  # Post-SCF
    assert cp2k._get_run_type(MockCalc(["pbe", "tddft"])) == "ENERGY"  # Post-SCF
    assert cp2k._get_run_type(MockCalc(["pbe"])) == "ENERGY"  # Default

    print("  ✓ Run type detection works correctly")


def test_extended_xyz_parsing():
    """Test extended XYZ trajectory parsing (CP2K 2026.1+)"""
    print("Testing extended XYZ trajectory parsing...")

    cp2k = CP2K()

    # Create a test extended XYZ file
    extended_xyz_content = """3
Properties=species:S:1:pos:R:3:forces:R:3 energy=-76.123456 pbc="F F F"
O   0.000000   0.000000   0.117000   0.001   0.002  -0.003
H   0.000000   0.757000  -0.469000   0.005   0.001   0.001
H   0.000000  -0.757000  -0.469000  -0.006  -0.003   0.002
3
Properties=species:S:1:pos:R:3:forces:R:3 energy=-76.234567 pbc="F F F"
O   0.000000   0.000000   0.120000   0.000   0.000   0.000
H   0.000000   0.760000  -0.470000   0.000   0.000   0.000
H   0.000000  -0.760000  -0.470000   0.000   0.000   0.000
"""

    # Write to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        f.write(extended_xyz_content)
        temp_path = f.name

    try:
        # Parse the trajectory
        coords = cp2k._parse_xyz_trajectory(temp_path)

        assert coords is not None, "Failed to parse extended XYZ"
        assert len(coords) == 3, f"Expected 3 atoms, got {len(coords)}"

        # Check last frame coordinates (second frame)
        assert abs(coords[0][2] - 0.120) < 0.001, "O z-coordinate mismatch"
        assert abs(coords[1][1] - 0.760) < 0.001, "H1 y-coordinate mismatch"
        assert abs(coords[2][1] - (-0.760)) < 0.001, "H2 y-coordinate mismatch"

        print("  ✓ Extended XYZ trajectory parsing works correctly")
    finally:
        os.unlink(temp_path)


def test_standard_xyz_parsing():
    """Test standard XYZ trajectory parsing (backward compatibility)"""
    print("Testing standard XYZ trajectory parsing...")

    cp2k = CP2K()

    # Create a standard XYZ file
    standard_xyz_content = """3
Water molecule frame 1
O   0.000000   0.000000   0.117000
H   0.000000   0.757000  -0.469000
H   0.000000  -0.757000  -0.469000
3
Water molecule frame 2
O   0.000000   0.000000   0.115000
H   0.000000   0.755000  -0.468000
H   0.000000  -0.755000  -0.468000
"""

    # Write to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        f.write(standard_xyz_content)
        temp_path = f.name

    try:
        # Parse the trajectory
        coords = cp2k._parse_xyz_trajectory(temp_path)

        assert coords is not None, "Failed to parse standard XYZ"
        assert len(coords) == 3, f"Expected 3 atoms, got {len(coords)}"

        # Check last frame coordinates (second frame)
        assert abs(coords[0][2] - 0.115) < 0.001, "O z-coordinate mismatch"
        assert abs(coords[1][1] - 0.755) < 0.001, "H1 y-coordinate mismatch"

        print("  ✓ Standard XYZ trajectory parsing works correctly")
    finally:
        os.unlink(temp_path)


def test_dispersion_detection():
    """Test dispersion correction detection"""
    print("Testing dispersion correction detection...")

    cp2k = CP2K()

    class MockKeyword:
        def __init__(self, kw):
            self._kw = kw
        def __str__(self):
            return self._kw

    class MockInput:
        def __init__(self, keywords):
            self.keywords = [MockKeyword(k) for k in keywords]

    class MockCalc:
        def __init__(self, keywords):
            self.input = MockInput(keywords)

    # Test dispersion detection
    assert cp2k._has_dispersion(MockCalc(["pbe", "d3"])) == True
    assert cp2k._has_dispersion(MockCalc(["pbe", "d4"])) == True
    assert cp2k._has_dispersion(MockCalc(["pbe", "disp"])) == True
    assert cp2k._has_dispersion(MockCalc(["pbe"])) == False

    print("  ✓ Dispersion correction detection works correctly")


def test_functional_extraction():
    """Test functional extraction from keywords"""
    print("Testing functional extraction...")

    cp2k = CP2K()

    class MockKeyword:
        def __init__(self, kw):
            self._kw = kw
        def __str__(self):
            return self._kw

    class MockInput:
        def __init__(self, keywords):
            self.keywords = [MockKeyword(k) for k in keywords]

    class MockCalc:
        def __init__(self, keywords):
            self.input = MockInput(keywords)

    # Test functional extraction
    assert cp2k._get_functional(MockCalc(["pbe", "opt"])) == "pbe"
    assert cp2k._get_functional(MockCalc(["b3lyp", "opt"])) == "b3lyp"
    assert cp2k._get_functional(MockCalc(["pbe0", "opt"])) == "pbe0"
    assert cp2k._get_functional(MockCalc(["blyp", "opt"])) == "blyp"
    assert cp2k._get_functional(MockCalc(["tpss", "opt"])) == "tpss"
    assert cp2k._get_functional(MockCalc(["opt"])) == "pbe"  # Default

    print("  ✓ Functional extraction works correctly")


def test_gth_potential_mapping():
    """Test GTH pseudopotential mapping"""
    print("Testing GTH pseudopotential mapping...")

    cp2k = CP2K()

    # Test common elements
    assert cp2k._cp2k_potential("H") == "GTH-PBE-q1"
    assert cp2k._cp2k_potential("C") == "GTH-PBE-q4"
    assert cp2k._cp2k_potential("N") == "GTH-PBE-q5"
    assert cp2k._cp2k_potential("O") == "GTH-PBE-q6"
    assert cp2k._cp2k_potential("S") == "GTH-PBE-q6"
    assert cp2k._cp2k_potential("Fe") == "GTH-PBE-q16"
    assert cp2k._cp2k_potential("Pd") == "GTH-PBE-q18"
    assert cp2k._cp2k_potential("Au") == "GTH-PBE-q11"

    # Unknown element should default to q4
    assert cp2k._cp2k_potential("Uue") == "GTH-PBE-q4"

    print("  ✓ GTH pseudopotential mapping works correctly")


def run_all_tests():
    """Run all unit tests"""
    print("=" * 60)
    print("CP2K 2026.1 Feature Tests for autodE Wrapper")
    print("=" * 60)
    print()

    tests = [
        test_aug_molopt_basis_mapping,
        test_sftddft_detection,
        test_tddft_detection,
        test_nstates_parsing,
        test_cg_optimizer_detection,
        test_run_type_detection,
        test_extended_xyz_parsing,
        test_standard_xyz_parsing,
        test_dispersion_detection,
        test_functional_extraction,
        test_gth_potential_mapping,
    ]

    passed = 0
    failed = 0

    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"  ✗ FAILED: {e}")
            failed += 1
        except Exception as e:
            print(f"  ✗ ERROR: {e}")
            failed += 1

    print()
    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)

    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
