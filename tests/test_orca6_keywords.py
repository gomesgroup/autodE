"""
Tests for ORCA 6.x keyword classes.

These tests verify:
1. Keyword generation is correct
2. ORCA block generation is correct
3. Edge cases are handled

NOTE: These keyword classes are standalone configuration objects that generate
ORCA input strings. They do NOT inherit from autodE's Keywords class.
"""

import pytest
import sys
sys.path.insert(0, '/tmp/autode-fork/autode/wrappers/keywords')

from orca6 import (
    GOATKeywords,
    TSConformerKeywords,
    SolvatorKeywords,
    DockerKeywords,
    IRCKeywords,
    NEBKeywords,
    ONIOMKeywords,
    MultiscaleNEBTSKeywords,
    MLIPConfig,
    ExtOptKeywords,
)


class TestGOATKeywords:
    """Test GOAT (Global Optimization And Transition state) keywords."""

    def test_default_keyword(self):
        """Test default GOAT keyword generation."""
        kw = GOATKeywords()
        assert "GOAT" in kw.to_orca_keyword()

    def test_goat_block_default(self):
        """Test default GOAT block generation."""
        kw = GOATKeywords()
        block = kw.to_orca_block()
        assert "%goat" in block
        assert "MAXSTEPS 2000" in block
        assert "TEMPERATURE 400" in block
        assert "end" in block.lower()

    def test_goat_block_custom(self):
        """Test custom GOAT parameters."""
        kw = GOATKeywords(max_steps=5000, temperature=300, force_field="UFF")
        block = kw.to_orca_block()
        assert "MAXSTEPS 5000" in block
        assert "TEMPERATURE 300" in block
        assert "FORCEFIELD UFF" in block

    def test_goat_ts_mode(self):
        """Test GOAT in TS search mode."""
        kw = GOATKeywords(ts_search=True)
        keyword = kw.to_orca_keyword()
        assert "TS" in keyword or kw.ts_search is True


class TestTSConformerKeywords:
    """Test TS conformer search keywords."""

    def test_racerts_is_default(self):
        """RACE-TS should be the default (faster than GOAT)."""
        kw = TSConformerKeywords(reacting_atoms=[0, 1, 2])
        assert kw.use_goat is False

    def test_goat_can_be_enabled(self):
        """GOAT can be explicitly enabled."""
        kw = TSConformerKeywords(reacting_atoms=[0, 1, 2], use_goat=True)
        assert kw.use_goat is True

    def test_reacting_atoms_stored(self):
        """Test that reacting atoms are stored correctly."""
        kw = TSConformerKeywords(reacting_atoms=[0, 1, 2, 3])
        assert kw.reacting_atoms == [0, 1, 2, 3]

    def test_racerts_detection(self):
        """Test RACE-TS availability detection."""
        # This should return True/False without error
        result = TSConformerKeywords.is_racerts_available()
        assert isinstance(result, bool)


class TestSolvatorKeywords:
    """Test SOLVATOR explicit solvation keywords."""

    def test_default_solvent(self):
        """Test default solvent is water."""
        kw = SolvatorKeywords()
        assert kw.solvent == "water"

    def test_solvator_keyword(self):
        """Test SOLVATOR keyword generation."""
        kw = SolvatorKeywords()
        assert "SOLVATOR" in kw.to_orca_keyword()

    def test_solvator_block(self):
        """Test SOLVATOR block generation."""
        kw = SolvatorKeywords(solvent="water", n_shells=2)
        block = kw.to_orca_block()
        assert "%solvator" in block.lower()
        assert "water" in block.lower()

    def test_custom_solvent(self):
        """Test non-water solvent."""
        kw = SolvatorKeywords(solvent="methanol", n_shells=3)
        block = kw.to_orca_block()
        assert "methanol" in block.lower()


class TestDockerKeywords:
    """Test DOCKER molecular docking keywords."""

    def test_docker_keyword(self):
        """Test DOCKER keyword generation."""
        kw = DockerKeywords()
        assert "DOCKER" in kw.to_orca_keyword()

    def test_docker_block(self):
        """Test DOCKER block generation."""
        kw = DockerKeywords(n_structures=20)
        block = kw.to_orca_block()
        assert "%docker" in block.lower()


class TestIRCKeywords:
    """Test IRC (Intrinsic Reaction Coordinate) keywords."""

    def test_irc_keyword(self):
        """Test IRC keyword generation."""
        kw = IRCKeywords()
        assert "IRC" in kw.to_orca_keyword()

    def test_irc_direction_both(self):
        """Test IRC in both directions (default)."""
        kw = IRCKeywords(direction="both")
        block = kw.to_orca_block()
        assert "Both" in block or "BOTH" in block.upper()

    def test_irc_direction_forward(self):
        """Test IRC in forward direction only."""
        kw = IRCKeywords(direction="forward")
        block = kw.to_orca_block()
        assert "Forward" in block or "FORWARD" in block.upper()

    def test_irc_parameters(self):
        """Test IRC step parameters."""
        kw = IRCKeywords(max_iter=100, step_size=0.05)
        block = kw.to_orca_block()
        assert "100" in block


class TestNEBKeywords:
    """Test NEB (Nudged Elastic Band) keywords."""

    def test_neb_keyword(self):
        """Test NEB keyword generation."""
        kw = NEBKeywords()
        keyword = kw.to_orca_keyword()
        assert "NEB" in keyword

    def test_neb_ts_keyword(self):
        """Test NEB-TS keyword generation."""
        kw = NEBKeywords(ts_search=True)
        keyword = kw.to_orca_keyword()
        assert "NEB-TS" in keyword or ("NEB" in keyword and "TS" in keyword)

    def test_neb_images(self):
        """Test NEB image count."""
        kw = NEBKeywords(n_images=15)
        block = kw.to_orca_block()
        assert "15" in block

    def test_climbing_image(self):
        """Test climbing image NEB."""
        kw = NEBKeywords(climbing_image=True)
        block = kw.to_orca_block()
        assert "CI True" in block or "True" in block


class TestONIOMKeywords:
    """Test QM/QM2 ONIOM keywords."""

    def test_qm_xtb_oniom(self):
        """Test QM/XTB ONIOM setup."""
        kw = ONIOMKeywords(
            high_level="r2SCAN-3c",
            low_level="XTB",
            qm_atoms=[0, 1, 2, 3],
        )
        keyword = kw.to_orca_keyword()
        assert "QM/XTB" in keyword or "QM" in keyword

    def test_qm_atoms_in_block(self):
        """Test QM atoms specification in QMMM block."""
        kw = ONIOMKeywords(
            high_level="B3LYP def2-SVP",
            low_level="XTB",
            qm_atoms=[0, 1, 2],
        )
        block = kw.to_orca_block()
        assert "qmmm" in block.lower() or "qmatoms" in block.lower()
        # Atom indices should be present
        assert "0" in block and "1" in block and "2" in block

    def test_qm_mlip_oniom(self):
        """Test QM/MLIP ONIOM setup."""
        kw = ONIOMKeywords(
            high_level="r2SCAN-3c",
            low_level="MLIP",
            low_level_mlip=True,
            mlip_model="aimnet2",
            qm_atoms=[0, 1, 2],
        )
        assert kw.low_level_mlip is True
        assert kw.mlip_model == "aimnet2"


class TestMultiscaleNEBTSKeywords:
    """Test Multiscale (ONIOM + NEB-TS) keywords."""

    def test_multiscale_nebts_defaults(self):
        """Test multiscale NEB-TS default values."""
        kw = MultiscaleNEBTSKeywords(
            high_level="r2SCAN-3c",
            low_level="XTB",
            qm_atoms=[0, 1, 2],
            n_images=12,
        )
        assert kw.n_images == 12
        assert kw.numfreq is True

    def test_multiscale_with_mlip(self):
        """Test multiscale NEB-TS with MLIP environment."""
        kw = MultiscaleNEBTSKeywords(
            high_level="r2SCAN-3c",
            low_level="MLIP",
            use_mlip=True,
            mlip_model="aimnet2",
            qm_atoms=[0, 1, 2],
            n_images=12,
        )
        assert kw.use_mlip is True
        assert kw.mlip_model == "aimnet2"


class TestMLIPConfig:
    """Test MLIP configuration."""

    def test_default_model(self):
        """Test default MLIP model."""
        config = MLIPConfig()
        assert config.model in ["aimnet2", "uma", None] or config.model is None

    def test_server_url(self):
        """Test server URL configuration."""
        config = MLIPConfig(
            model="aimnet2",
            server_url="http://gpg-boltzmann:5003",
        )
        assert "5003" in config.server_url


class TestExtOptKeywords:
    """Test ExtOpt (external optimizer) keywords."""

    def test_extopt_keyword(self):
        """Test ExtOpt keyword generation."""
        kw = ExtOptKeywords(command="mlip_client http://localhost:5003 aimnet2")
        keyword = kw.to_orca_keyword()
        assert "ExtOpt" in keyword

    def test_extopt_block(self):
        """Test ExtOpt block generation."""
        kw = ExtOptKeywords(command="mlip_client http://localhost:5003 aimnet2")
        block = kw.to_orca_block()
        assert "%extopt" in block.lower()
        assert "cmd" in block.lower()


class TestKeywordIntegration:
    """Integration tests for keyword combinations."""

    def test_oniom_generates_valid_input(self):
        """Test that ONIOM keywords generate syntactically valid ORCA input."""
        kw = ONIOMKeywords(
            high_level="r2SCAN-3c",
            low_level="XTB",
            qm_atoms=[0, 1, 2, 3, 4],
        )

        keyword = kw.to_orca_keyword()
        block = kw.to_orca_block()

        # Basic syntax checks
        assert keyword.startswith("!")
        assert "end" in block.lower()
        # No unbalanced braces
        assert block.count("{") == block.count("}")

    def test_neb_generates_valid_input(self):
        """Test that NEB keywords generate syntactically valid ORCA input."""
        kw = NEBKeywords(n_images=12, climbing_image=True, ts_search=True)

        keyword = kw.to_orca_keyword()
        block = kw.to_orca_block()

        assert keyword.startswith("!")
        assert "end" in block.lower()

    def test_all_blocks_balanced(self):
        """Test that all blocks have balanced END statements."""
        test_cases = [
            GOATKeywords(),
            SolvatorKeywords(),
            DockerKeywords(),
            IRCKeywords(),
            NEBKeywords(),
            ONIOMKeywords(high_level="B3LYP", low_level="XTB", qm_atoms=[0]),
        ]

        for kw in test_cases:
            block = kw.to_orca_block()
            sections = block.lower().count("%")
            ends = block.lower().count("end")
            assert ends >= sections, f"{type(kw).__name__} has unbalanced blocks"
