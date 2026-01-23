"""
Tests for ORCA 6.x wrapper functions.

These tests verify:
1. Input generation is correct
2. Output parsing works
3. API functions work as expected
"""

import pytest
import os
import sys
import tempfile
from unittest.mock import patch, MagicMock

sys.path.insert(0, '/tmp/autode-fork/autode/wrappers')
sys.path.insert(0, '/tmp/autode-fork/autode/wrappers/keywords')


class TestSolvatorWrapper:
    """Test SOLVATOR wrapper functions."""

    def test_generate_solvator_input(self):
        """Test SOLVATOR input generation."""
        from solvator import generate_solvator_input
        from orca6 import SolvatorKeywords

        kw = SolvatorKeywords(solvent="water", n_shells=2)
        coords = [
            ("O", 0.0, 0.0, 0.0),
            ("H", 0.0, 0.8, 0.6),
            ("H", 0.0, -0.8, 0.6),
        ]

        input_str = generate_solvator_input(
            keywords=kw,
            charge=0,
            multiplicity=1,
            coordinates=coords,
            method="B3LYP",
            basis="def2-SVP",
            n_cores=4,
        )

        # Check basic structure
        assert "SOLVATOR" in input_str or "solvator" in input_str.lower()
        assert "B3LYP" in input_str
        assert "def2-SVP" in input_str
        assert "*xyz 0 1" in input_str
        assert "nprocs 4" in input_str

    def test_parse_solvator_output(self):
        """Test SOLVATOR output parsing."""
        from solvator import parse_solvator_output

        # Create mock output file
        mock_output = """
ORCA 6.1.1

Running SOLVATOR calculation...
30 solvent molecules added

FINAL SINGLE POINT ENERGY      -123.456789

ORCA TERMINATED NORMALLY
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.out', delete=False) as f:
            f.write(mock_output)
            temp_path = f.name

        try:
            result = parse_solvator_output(temp_path)
            assert result["terminated_normally"] is True
            assert result["n_solvent_molecules"] == 30
            assert abs(result["solvation_energy"] - (-123.456789)) < 1e-6
        finally:
            os.unlink(temp_path)


class TestDockerWrapper:
    """Test DOCKER wrapper functions."""

    def test_generate_docker_input(self):
        """Test DOCKER input generation."""
        from docker import generate_docker_input
        from orca6 import DockerKeywords

        kw = DockerKeywords(n_structures=10)
        mol1_coords = [("Cl", -2.0, 0.0, 0.0)]
        mol2_coords = [
            ("C", 0.0, 0.0, 0.0),
            ("H", 1.0, 0.0, 0.0),
            ("Br", 2.0, 0.0, 0.0),
        ]

        input_str = generate_docker_input(
            keywords=kw,
            charge=-1,
            multiplicity=1,
            mol1_coords=mol1_coords,
            mol2_coords=mol2_coords,
            method="B3LYP",
            basis="def2-SVP",
            n_cores=8,
        )

        assert "DOCKER" in input_str.upper()
        assert "*xyz -1 1" in input_str
        assert "Cl" in input_str
        assert "Br" in input_str

    def test_parse_docker_output(self):
        """Test DOCKER output parsing."""
        from docker import parse_docker_output

        mock_output = """
ORCA 6.1.1

Running DOCKER calculation...
10 structures generated

FINAL SINGLE POINT ENERGY      -456.789012

ORCA TERMINATED NORMALLY
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.out', delete=False) as f:
            f.write(mock_output)
            temp_path = f.name

        try:
            result = parse_docker_output(temp_path)
            assert result["terminated_normally"] is True
            assert result["n_structures_found"] == 10
        finally:
            os.unlink(temp_path)


class TestIRCWrapper:
    """Test IRC wrapper functions."""

    def test_generate_irc_input(self):
        """Test IRC input generation."""
        from irc import generate_irc_input
        from orca6 import IRCKeywords

        kw = IRCKeywords(direction="both", max_iter=50)
        coords = [
            ("C", 0.0, 0.0, 0.0),
            ("H", 1.0, 0.0, 0.0),
            ("H", -1.0, 0.0, 0.0),
        ]

        input_str = generate_irc_input(
            keywords=kw,
            charge=0,
            multiplicity=1,
            coordinates=coords,
            method="B3LYP",
            basis="def2-SVP",
            n_cores=8,
        )

        assert "IRC" in input_str.upper()
        assert "*xyz 0 1" in input_str

    def test_parse_irc_output(self):
        """Test IRC output parsing."""
        from irc import parse_irc_output

        mock_output = """
ORCA 6.1.1

Starting IRC calculation...

FINAL SINGLE POINT ENERGY      -100.123456

IRC FORWARD
Step 1: Reaction coordinate = 0.1
FINAL SINGLE POINT ENERGY      -100.123000

IRC BACKWARD
Step 1: Reaction coordinate = -0.1
FINAL SINGLE POINT ENERGY      -100.124000

IRC FINISHED

ORCA TERMINATED NORMALLY
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.out', delete=False) as f:
            f.write(mock_output)
            temp_path = f.name

        try:
            result = parse_irc_output(temp_path)
            assert result.ts_energy is not None
            assert result.connects_to_products or len(result.forward_path.energies) > 0
            assert result.connects_to_reactants or len(result.backward_path.energies) > 0
        finally:
            os.unlink(temp_path)

    def test_calculate_activation_energies(self):
        """Test activation energy calculation."""
        from irc import calculate_activation_energies, IRCResult, IRCPath

        result = IRCResult(
            ts_structure=[("H", 0, 0, 0)],
            ts_energy=-100.0,
            forward_path=IRCPath(energies=[-100.5]),
            backward_path=IRCPath(energies=[-100.8]),
        )

        fwd_barrier, rev_barrier = calculate_activation_energies(result)

        # TS - reactant (backward endpoint)
        assert abs(fwd_barrier - 0.8) < 1e-6
        # TS - product (forward endpoint)
        assert abs(rev_barrier - 0.5) < 1e-6


class TestMLIPWrapper:
    """Test MLIP external wrapper functions."""

    def test_check_mlip_server_timeout(self):
        """Test MLIP server check handles timeout gracefully."""
        from mlip_external import check_mlip_server

        # Non-existent server should return False, not raise
        result = check_mlip_server("http://nonexistent.server:9999")
        assert result is False

    def test_find_best_mlip_server_no_server(self):
        """Test find_best_mlip_server returns None when no servers available."""
        from mlip_external import find_best_mlip_server

        with patch('mlip_external.check_mlip_server', return_value=False):
            result = find_best_mlip_server()
            assert result is None

    def test_create_extopt_script(self):
        """Test ExtOpt script creation."""
        from mlip_external import create_extopt_script

        with tempfile.TemporaryDirectory() as tmpdir:
            script_path = os.path.join(tmpdir, "test_extopt.sh")
            result = create_extopt_script(
                model="aimnet2",
                server_url="http://localhost:5003",
                output_path=script_path,
            )

            assert os.path.exists(result)
            with open(result, 'r') as f:
                content = f.read()
            assert "#!/bin/bash" in content
            assert "aimnet2" in content
            assert "localhost:5003" in content

    def test_generate_qm_mlip_oniom_input(self):
        """Test QM/MLIP ONIOM input generation."""
        from mlip_external import generate_qm_mlip_oniom_input

        coords = [
            ("C", 0.0, 0.0, 0.0),
            ("H", 1.0, 0.0, 0.0),
            ("H", -1.0, 0.0, 0.0),
            ("H", 0.0, 1.0, 0.0),
            ("H", 0.0, -1.0, 0.0),
        ]

        input_str = generate_qm_mlip_oniom_input(
            high_level_method="r2SCAN-3c",
            mlip_model="aimnet2",
            qm_atoms=[0, 1, 2],
            coordinates=coords,
            charge=0,
            multiplicity=1,
            server_url="http://localhost:5003",
            n_cores=8,
        )

        # Should have QM/QM2 keyword
        assert "QM/QM2" in input_str or "ONIOM" in input_str.upper()
        # Should reference ExtOpt
        assert "extopt" in input_str.lower() or "ExtOpt" in input_str
        # Should have QM atoms
        assert "0" in input_str and "1" in input_str and "2" in input_str
        # Should have coordinates
        assert "*xyz 0 1" in input_str


class TestMLIPAcceleratedNEB:
    """Test MLIP-accelerated NEB class."""

    def test_initialization(self):
        """Test MLIPAcceleratedNEB initialization."""
        from mlip_external import MLIPAcceleratedNEB

        # Create mock molecules
        mock_reactant = MagicMock()
        mock_product = MagicMock()

        neb = MLIPAcceleratedNEB(
            reactant=mock_reactant,
            product=mock_product,
            mlip_model="aimnet2",
            dft_method="r2SCAN-3c",
            n_images=15,
        )

        assert neb.mlip_model == "aimnet2"
        assert neb.dft_method == "r2SCAN-3c"
        assert neb.n_images == 15
        assert neb.ts_guess is None  # Not computed yet


class TestInputSyntax:
    """Test that generated inputs have valid ORCA syntax."""

    def test_no_duplicate_keywords(self):
        """Test that keywords don't have duplicates."""
        from orca6 import ONIOMKeywords

        kw = ONIOMKeywords(
            high_level="r2SCAN-3c",
            low_level="XTB",
            qm_atoms=[0, 1],
        )

        keyword = kw.to_orca_keyword()
        # Check no duplicate words
        words = keyword.replace("!", "").split()
        assert len(words) == len(set(words)), f"Duplicate keywords in: {keyword}"

    def test_balanced_blocks(self):
        """Test that blocks have balanced END statements."""
        from orca6 import (
            GOATKeywords,
            SolvatorKeywords,
            DockerKeywords,
            IRCKeywords,
            NEBKeywords,
            ONIOMKeywords,
        )

        for KeywordClass in [GOATKeywords, SolvatorKeywords, DockerKeywords,
                              IRCKeywords, NEBKeywords, ONIOMKeywords]:
            if KeywordClass == ONIOMKeywords:
                kw = KeywordClass(high_level="B3LYP", low_level="XTB", qm_atoms=[0])
            else:
                kw = KeywordClass()

            block = kw.to_orca_block()
            # Count % sections vs end statements
            sections = block.lower().count("%")
            ends = block.lower().count("end")
            assert ends >= sections, f"{KeywordClass.__name__} has unbalanced blocks"
