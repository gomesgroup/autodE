import autode as ade

# Create an aspirin molecule
m = ade.Molecule(smiles='CC(=O)Oc1ccccc1C(=O)O')

print(m.fragments)

m.fragments.fragment_constraints = {
    1: {
        "atom_idxs": [0, 1, 2, 3, 4],
        "strategy": "fix"
    }
}

print(m.fragments)

# Test ORCA implementation
m.optimise(method=ade.methods.ORCA(), n_cores=32)
