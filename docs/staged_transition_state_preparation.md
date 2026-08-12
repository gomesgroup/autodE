# Staged transition-state preparation

`autode.transition_states.StagedTSPreparation` provides a method-agnostic
four-stage preparation for a transition-state guess:

1. hold the active atoms fixed and relax the surrounding structure;
2. release the active atoms while holding target active-bond distances;
3. hold spectator atoms and run an active-region transition-state search;
4. remove all constraints and run a final full transition-state search.

The caller supplies one exact autoDE `Method` instance.  The class never
resolves `Config.hcode` or `Config.lcode`.  A `calculation_runner` hook receives
every autoDE `Calculation` before execution, allowing a project to inject and
audit engine-specific controls, preserve provenance, or route calculations to
a remote potential.

```python
from autode.transition_states import StagedTSPreparation

preparation = StagedTSPreparation(
    method=explicit_method,
    active_atom_indices=(carbon, nitrogen, leaving_group),
    active_bond_distances={(carbon, nitrogen): 2.0, (carbon, leaving_group): 2.3},
    calculation_runner=audited_runner,
    n_cores=4,
)
result = preparation.run(ts_guess)
transition_state = result.transition_state
```

For complete autoDE transition-state location, exact high- and low-level
methods can also be supplied directly:

```python
reaction.locate_transition_state(
    hmethod=explicit_method,
    lmethod=explicit_method,
)
```

Exact `Method` objects are also valid values for `Config.hcode` and
`Config.lcode`.  Explicit arguments are preferred for reproducible workflows
because they are propagated through template searches, AdaptivePath, OptTS,
mode validation, endpoint relaxation, and TS conformer refinement.

`MLIPAcceleratedNEB` also retains the current peak as a labelled TS guess when
a climbing-image band stops at its iteration budget.  The result remains
`nonconverged`; downstream code can refine the seed without confusing it with a
converged path.

For external-I/O CINEB calculations, `run_mlip_neb` accepts independent
`image_workers` and `calculation_cores` values and rejects requests whose
product exceeds the declared total core allocation.  This permits, for
example, four concurrent one-core ProgExt gradients inside a four-core job.
