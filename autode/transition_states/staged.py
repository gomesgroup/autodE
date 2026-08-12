"""Reusable constrained preparation for transition-state optimisation.

The implementation is deliberately method-agnostic.  Projects inject one
exact :class:`Method` instance and may wrap every :class:`Calculation` through
``calculation_runner`` to audit or modify generated input at the executor
boundary.  No configured high- or low-level method is resolved here.
"""

from dataclasses import dataclass
from typing import Any, Callable, Mapping, Optional, Sequence, Tuple, TYPE_CHECKING

from autode.calculations import Calculation
from autode.config import Config
from autode.constraints import DistanceConstraints
from autode.transition_states.transition_state import TransitionState

if TYPE_CHECKING:
    from autode.transition_states.ts_guess import TSguess
    from autode.wrappers.keywords import Keywords
    from autode.wrappers.methods import Method


StageRunner = Callable[[Calculation], Any]


@dataclass(frozen=True)
class StagedTSKeywords:
    """Calculation keywords for the four preparation stages."""

    environment: Optional["Keywords"] = None
    active_bonds: Optional["Keywords"] = None
    active_region_ts: Optional["Keywords"] = None
    final_ts: Optional["Keywords"] = None


@dataclass(frozen=True)
class StagedTSReceipt:
    """One completed calculation in a staged TS preparation."""

    stage: str
    calculation_name: str
    method_name: str
    n_cores: int
    fixed_atom_indices: Tuple[int, ...]
    distance_constraints: Tuple[Tuple[int, int, float], ...]
    energy_hartree: Optional[float]
    converged: Optional[bool]


@dataclass(frozen=True)
class StagedTSResult:
    """Final TS candidate and immutable stage receipts."""

    transition_state: TransitionState
    receipts: Tuple[StagedTSReceipt, ...]


class StagedTSPreparationError(RuntimeError):
    """Raised when a required staged optimisation does not converge."""


class StagedTSPreparation:
    """Prepare a TS through environment, bond, active, and full stages.

    The active atom set usually contains atoms participating in forming or
    breaking bonds.  Stage 1 fixes those atoms while the environment relaxes;
    stage 2 releases them but fixes active-bond distances; stage 3 fixes all
    spectator atoms for an active-region OptTS; and stage 4 removes every
    constraint for a full OptTS.
    """

    _STAGE_NAMES = (
        "environment_relax",
        "active_bond_relax",
        "active_region_optts",
        "full_optts",
    )

    def __init__(
        self,
        *,
        method: "Method",
        active_atom_indices: Sequence[int],
        active_bond_distances: Optional[Mapping[Tuple[int, int], Any]] = None,
        keywords: Optional[StagedTSKeywords] = None,
        n_cores: Optional[int] = None,
        calculation_runner: Optional[StageRunner] = None,
        name_prefix: str = "staged_ts",
        required_converged_stages: Sequence[str] = ("full_optts",),
    ) -> None:
        if method is None:
            raise ValueError("Staged TS preparation requires an exact method")
        self.method = method
        self.active_atom_indices = tuple(
            sorted({int(index) for index in active_atom_indices})
        )
        if not self.active_atom_indices:
            raise ValueError("At least one active atom is required")
        if any(index < 0 for index in self.active_atom_indices):
            raise ValueError("Active atom indices must be non-negative")

        self.active_bond_distances = (
            None
            if active_bond_distances is None
            else DistanceConstraints(dict(active_bond_distances))
        )
        self.keywords = keywords or StagedTSKeywords()
        self.n_cores = Config.n_cores if n_cores is None else int(n_cores)
        if self.n_cores < 1:
            raise ValueError("Staged TS preparation requires at least one core")
        self.calculation_runner = calculation_runner or (lambda calc: calc.run())
        self.name_prefix = str(name_prefix)
        required = tuple(dict.fromkeys(map(str, required_converged_stages)))
        unknown = set(required).difference(self._STAGE_NAMES)
        if unknown:
            raise ValueError(
                "Unknown required staged TS operations: "
                + ", ".join(sorted(unknown))
            )
        self.required_converged_stages = frozenset(required)

    def _keywords_for(self, stage: str):
        if stage == "environment_relax":
            return (
                self.keywords.environment
                or self.method.keywords.opt
                or self.method.keywords.low_opt
            )
        if stage == "active_bond_relax":
            return self.keywords.active_bonds or self.method.keywords.opt
        if stage == "active_region_optts":
            return self.keywords.active_region_ts or self.method.keywords.opt_ts
        if stage == "full_optts":
            return self.keywords.final_ts or self.method.keywords.opt_ts
        raise ValueError(f"Unknown staged TS operation {stage}")

    def _run_stage(
        self,
        molecule,
        stage: str,
        fixed_atoms: Sequence[int] = (),
        distances: Optional[Mapping[Tuple[int, int], Any]] = None,
    ) -> StagedTSReceipt:
        keywords = self._keywords_for(stage)
        if keywords is None:
            raise ValueError(
                f"Method {self.method.name} has no keywords for {stage}"
            )

        molecule.constraints.cartesian = list(fixed_atoms)
        molecule.constraints.distance = (
            None if distances is None else dict(distances)
        )
        calc = Calculation(
            name=f"{self.name_prefix}_{stage}",
            molecule=molecule,
            method=self.method,
            keywords=keywords,
            n_cores=self.n_cores,
        )
        try:
            self.calculation_runner(calc)
        finally:
            molecule.constraints.cartesian = None
            molecule.constraints.distance = None

        constrained_distances = ()
        if distances is not None:
            constrained_distances = tuple(
                (int(pair[0]), int(pair[1]), float(distance))
                for pair, distance in sorted(
                    dict(distances).items(), key=lambda item: tuple(item[0])
                )
            )
        energy = molecule.energy
        try:
            converged = bool(calc.optimiser.converged)
        except (AttributeError, RuntimeError, ValueError):
            converged = None
        return StagedTSReceipt(
            stage=stage,
            calculation_name=calc.name,
            method_name=self.method.name,
            n_cores=self.n_cores,
            fixed_atom_indices=tuple(fixed_atoms),
            distance_constraints=constrained_distances,
            energy_hartree=None if energy is None else float(energy),
            converged=converged,
        )

    def run(self, ts_guess: "TSguess") -> StagedTSResult:
        """Run all four stages and return the final unconstrained TS.

        Conditioning stages may intentionally stop at bounded iteration limits.
        By default only the final unconstrained OptTS must report convergence;
        callers can require any additional stages through
        ``required_converged_stages``.
        """

        if max(self.active_atom_indices) >= ts_guess.n_atoms:
            raise ValueError("Active atom index exceeds TS atom count")

        working = ts_guess.copy()
        receipts = []

        def append_stage(receipt: StagedTSReceipt) -> None:
            receipts.append(receipt)
            if (
                receipt.stage in self.required_converged_stages
                and receipt.converged is not True
            ):
                raise StagedTSPreparationError(
                    f"Required staged TS operation {receipt.stage} did not "
                    f"converge (reported {receipt.converged})"
                )

        append_stage(
            self._run_stage(
                working, "environment_relax", fixed_atoms=self.active_atom_indices
            )
        )

        distances = self.active_bond_distances
        if distances is None:
            if working.bond_rearrangement is None:
                raise ValueError(
                    "Active bond distances or a bond rearrangement are required"
                )
            distances = DistanceConstraints(
                {
                    tuple(pair): working.distance(*pair)
                    for pair in working.bond_rearrangement.all
                }
            )
        append_stage(
            self._run_stage(
                working,
                "active_bond_relax",
                distances=distances,
            )
        )

        transition_state = TransitionState(
            working, bond_rearr=working.bond_rearrangement
        )
        active = set(self.active_atom_indices)
        spectators = tuple(
            index for index in range(transition_state.n_atoms) if index not in active
        )
        append_stage(
            self._run_stage(
                transition_state,
                "active_region_optts",
                fixed_atoms=spectators,
            )
        )
        append_stage(self._run_stage(transition_state, "full_optts"))
        return StagedTSResult(
            transition_state=transition_state,
            receipts=tuple(receipts),
        )


__all__ = [
    "StagedTSKeywords",
    "StagedTSPreparation",
    "StagedTSPreparationError",
    "StagedTSReceipt",
    "StagedTSResult",
]
