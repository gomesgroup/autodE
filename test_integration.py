#!/usr/bin/env python3
"""
End-to-end integration test for ORCA 6.x keywords with autodE.

This tests that:
1. ORCA6Keywords can be used in a Calculation
2. Input files are generated correctly with both keyword lines AND blocks
3. All ORCA 6.x features integrate properly
"""
import sys
import os
sys.path.insert(0, '/tmp/autode-fork')

from autode import Molecule, Calculation
from autode.wrappers.ORCA import orca
from autode.wrappers.keywords import OptKeywords, SinglePointKeywords
from autode.wrappers.keywords.orca6 import (
    GOATKeywords,
    NEBKeywords,
    ONIOMKeywords,
    IRCKeywords,
    SolvatorKeywords,
    DockerKeywords,
    MultiscaleNEBTSKeywords,
    TSConformerKeywords,
)

def test_goat_integration():
    """Test GOAT optimization integration."""
    print("Testing GOAT optimization integration...")

    # Create a simple water molecule
    mol = Molecule(smiles='O')

    # Create GOAT keywords
    goat_kw = GOATKeywords(
        goat_type="OPT",
        method="r2SCAN-3c",
        max_steps=1000,
        temperature=300.0,
        force_field="GFN-FF",
        ewin=10.0,
    )

    # Create calculation with Keywords container
    calc = Calculation(
        name='goat_test',
        molecule=mol,
        method=orca,
        keywords=OptKeywords([goat_kw]),
        n_cores=4
    )

    # Generate input file
    calc.generate_input()

    # Read and verify input file
    with open('goat_test_orca.inp', 'r') as f:
        content = f.read()

    print("Generated input file:")
    print(content)
    print()

    # Verify keyword line
    assert 'GOAT_OPT r2SCAN-3c' in content, "GOAT keyword not in input file"

    # Verify block section
    assert '%goat' in content.lower(), "%goat block not found"
    assert 'MAXSTEPS 1000' in content, "MAXSTEPS not set correctly"
    assert 'TEMPERATURE 300.0' in content, "TEMPERATURE not set correctly"
    assert 'FORCEFIELD GFN-FF' in content, "FORCEFIELD not set correctly"
    assert 'EWIN 10.0' in content, "EWIN not set correctly"
    assert 'end' in content, "Block not closed with 'end'"

    print("✅ GOAT integration test PASSED")
    return True

def test_neb_integration():
    """Test NEB-TS integration."""
    print("\nTesting NEB-TS integration...")

    mol = Molecule(smiles='O')

    neb_kw = NEBKeywords(
        n_images=12,
        climbing_image=True,
        ts_search=True,
        spring_constant=0.05,
    )

    calc = Calculation(
        name='neb_test',
        molecule=mol,
        method=orca,
        keywords=OptKeywords([neb_kw]),
        n_cores=4
    )

    calc.generate_input()

    with open('neb_test_orca.inp', 'r') as f:
        content = f.read()

    print("Generated input file:")
    print(content)
    print()

    assert 'NEB-TS' in content, "NEB-TS keyword not in input file"
    assert '%neb' in content.lower(), "%neb block not found"
    assert 'NImages 12' in content, "NImages not set correctly"
    assert 'end' in content, "Block not closed"

    print("✅ NEB-TS integration test PASSED")
    return True

def test_oniom_integration():
    """Test QM/QM2 ONIOM integration."""
    print("\nTesting QM/QM2 ONIOM integration...")

    mol = Molecule(smiles='CCO')  # Ethanol

    oniom_kw = ONIOMKeywords(
        high_level="r2SCAN-3c",
        low_level="XTB",
        qm_atoms=[0, 1, 2],  # First 3 atoms in QM region
    )

    calc = Calculation(
        name='oniom_test',
        molecule=mol,
        method=orca,
        keywords=OptKeywords([oniom_kw]),
        n_cores=4
    )

    calc.generate_input()

    with open('oniom_test_orca.inp', 'r') as f:
        content = f.read()

    print("Generated input file:")
    print(content)
    print()

    assert 'QM' in content or 'ONIOM' in content.upper(), "ONIOM keyword not in input file"
    assert '%qmmm' in content.lower(), "%qmmm block not found"
    assert 'end' in content, "Block not closed"

    print("✅ ONIOM integration test PASSED")
    return True

def test_tsconformer_default():
    """Test that TSConformerKeywords defaults to RACE-TS (not GOAT)."""
    print("\nTesting TSConformerKeywords defaults to RACE-TS...")

    ts_kw = TSConformerKeywords(reacting_atoms=[0, 1, 2])

    assert ts_kw.use_goat is False, "TSConformerKeywords should default to use_goat=False (RACE-TS)"
    print("✅ TSConformerKeywords defaults to RACE-TS (use_goat=False)")

    # Test explicit GOAT mode
    ts_kw_goat = TSConformerKeywords(reacting_atoms=[0, 1, 2], use_goat=True)
    assert ts_kw_goat.use_goat is True, "TSConformerKeywords with use_goat=True should enable GOAT"
    print("✅ TSConformerKeywords can explicitly enable GOAT (use_goat=True)")

    return True

def test_solvator_integration():
    """Test SOLVATOR explicit solvation integration."""
    print("\nTesting SOLVATOR integration...")

    mol = Molecule(smiles='O')

    solvator_kw = SolvatorKeywords(
        solvent="water",
        n_shells=2,
    )

    calc = Calculation(
        name='solvator_test',
        molecule=mol,
        method=orca,
        keywords=OptKeywords([solvator_kw]),
        n_cores=4
    )

    calc.generate_input()

    with open('solvator_test_orca.inp', 'r') as f:
        content = f.read()

    print("Generated input file (first 500 chars):")
    print(content[:500])
    print("...\n")

    assert 'SOLVATOR' in content.upper(), "SOLVATOR keyword not in input file"
    assert '%solvator' in content.lower(), "%solvator block not found"
    assert 'end' in content, "Block not closed"

    print("✅ SOLVATOR integration test PASSED")
    return True

def run_all_tests():
    """Run all integration tests."""
    print("="*70)
    print("ORCA 6.x autodE Integration Tests")
    print("="*70)

    tests = [
        test_goat_integration,
        test_neb_integration,
        test_oniom_integration,
        test_tsconformer_default,
        test_solvator_integration,
    ]

    results = []
    for test in tests:
        try:
            result = test()
            results.append((test.__name__, result))
        except Exception as e:
            print(f"❌ {test.__name__} FAILED with error:")
            print(f"   {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
            results.append((test.__name__, False))

    print("\n" + "="*70)
    print("Test Results:")
    print("="*70)
    for name, passed in results:
        status = "✅ PASSED" if passed else "❌ FAILED"
        print(f"{name}: {status}")

    all_passed = all(r[1] for r in results)
    print("="*70)
    if all_passed:
        print("🎉 ALL TESTS PASSED! ORCA 6.x integration is fully functional.")
    else:
        print("⚠️  SOME TESTS FAILED. Review errors above.")
    print("="*70)

    return all_passed

if __name__ == '__main__':
    success = run_all_tests()
    sys.exit(0 if success else 1)
