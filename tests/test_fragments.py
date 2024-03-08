import autode as ade

# Create an aspirin molecule
m = ade.Molecule(smiles='CC(=O)Oc1ccccc1C(=O)O')

print(m.fragments)

# Use atom indices to define a fragment
m.fragments.fragment_constraints = {
    1: {
        "atom_idxs": [0, 1, 2, 3, 4],
        "strategy": "fix"
    }
}
print(m.fragments)

# Use start and end atoms to define a fragment
m.fragments.fragment_constraints = {
    1: {
        "start": 0,
        "end": 4,
        "strategy": "fix"
    },
    2: {
        "start": 5,
        "end": 9,
        "strategy": "relax"
    }
}
print(m.fragments)

# Test ORCA implementation
m.optimise(method=ade.methods.ORCA(), n_cores=32)
