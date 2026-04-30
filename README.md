[![Build Status](https://github.com/duartegroup/autodE/actions/workflows/pytest.yml/badge.svg)](https://github.com/duartegroup/autodE/actions) [![codecov](https://codecov.io/gh/duartegroup/autodE/branch/master/graph/badge.svg)](https://codecov.io/gh/duartegroup/autodE/branch/master) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![GitHub CodeQL](https://github.com/duartegroup/autodE/actions/workflows/codeql.yml/badge.svg)](https://github.com/duartegroup/autodE/actions/workflows/codeql.yml) [![Conda Recipe](https://img.shields.io/badge/recipe-autode-green.svg)](https://anaconda.org/conda-forge/autode) [![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/autode.svg)](https://anaconda.org/conda-forge/autode)

![alt text](autode/common/llogo.png)
***

## Introduction

**autodE** is a Python module initially designed for the automated calculation of reaction profiles from SMILES strings of
reactant(s) and product(s). Current features include: transition state location, conformer searching, atom mapping,
Python wrappers for a range of electronic structure theory codes, SMILES parsing, association complex generation, and
reaction profile generation.

## Why ORCA 6.x is the Recommended Backend

Starting with autodE v1.4.5 (gomesgroup fork), **ORCA 6.x is the recommended high-level QM backend** for automated reaction discovery. Here's why:

| Feature | Benefit for autodE | Alternative |
|---------|-------------------|-------------|
| **GOAT** (Global Optimizer with Adaptive Thermostat) | Automated conformer search with metadynamics - finds global minima without manual intervention | CREST (requires separate setup) |
| **SOLVATOR** | Automatic explicit solvation shell generation - critical for accurate solution-phase barriers | Manual placement |
| **DOCKER** | Built-in molecular docking for bimolecular reactions - automatic association complex generation | External docking tools |
| **NEB-TS** | Improved nudged elastic band for complex reaction paths | Standard NEB |
| **IRC** | Enhanced intrinsic reaction coordinate with bidirectional tracing | Manual IRC setup |
| **ExtOpt** | External optimizer interface - enables MLIP-accelerated geometry optimization | Pure DFT (slower) |
| **Free & Open** | No license fees, unlimited parallel jobs | Commercial alternatives |

### ORCA 6.x Feature Highlights

**GOAT (Metadynamics-Based Conformer Search)**
```python
import autode as ade
from autode.wrappers.keywords import GOATKeywords

# Automatic global minimum search
mol = ade.Molecule(smiles='CC(C)CC(C)C')  # 2,4-dimethylpentane
mol.optimise(method=ade.methods.ORCA(), keywords=GOATKeywords())
# GOAT automatically explores conformational space
```

**SOLVATOR (Explicit Solvation)**
```python
from autode.wrappers.keywords import SolvatorKeywords

# Add explicit water molecules around solute
keywords = SolvatorKeywords(solvent='water', n_solvent=50, shell_radius=5.0)
solvated_mol = mol.optimise(method=ade.methods.ORCA(), keywords=keywords)
```

**MLIP-Accelerated Optimization**
```python
from autode.wrappers import MLIPConfig

# Use machine learning potentials for fast pre-optimization
config = MLIPConfig(model='aimnet2', server_url='http://localhost:5003')
mol.optimise(method=ade.methods.ORCA(), mlip_config=config)
# 10-100x faster initial optimization before DFT refinement
```

**DOCKER (Bimolecular Complex Generation)**
```python
from autode.wrappers import dock_molecules

# Automatic docking for reaction complex setup
reactant1 = ade.Molecule(smiles='CCBr')
reactant2 = ade.Molecule(smiles='[OH-]')
complex = dock_molecules(reactant1, reactant2, method=ade.methods.ORCA())
# Finds optimal pre-reactive geometry for SN2 reaction
```

### g-xTB Driver

The fork ships a first-class `GXTB` low-level method (`autode/wrappers/GXTB.py`). It is registered in `autode.methods` so it can be selected by name (`Config.lcode = "gxtb"`) or imported directly:

```python
import autode as ade
from autode.methods import GXTB

# Single-point or optimization
mol = ade.Molecule(smiles="CC(=O)Oc1ccccc1C(=O)O", name="aspirin")
mol.optimise(method=GXTB())
print(mol.energy)

# As the low-level method in a reaction profile
ade.Config.lcode = "gxtb"
ade.Config.hcode = "orca"
rxn = ade.Reaction("CC(=O)O.CO>>CC(=O)OC.O", name="esterification")
rxn.calculate_reaction_profile()
```

`GXTB` is a thin subclass of `XTB` that:

- points at the g-xTB binary (override via the `GXTB_PATH` env var; defaults to `/mnt/beegfs/software/g-xtb-2.0.0/x86_64/bin/xtb` on the gpg cluster),
- forces `gfn_version = None` (g-xTB has its own parameter set),
- prepends `--gxtb` to the runtime flag list, and
- preserves the implicit `--opt` / `--grad` flags that the parent class would otherwise suppress when a non-empty keyword list is present.

Verified bit-identical SP energies vs. direct `xtb --gxtb` invocation; clean opt convergence on both x86_64 (native) and ARM64 (QEMU).

### Version Detection Fix

This fork fixes a critical bug where ORCA 6.x was incorrectly detected as ORCA 4.x:

```python
# Before (broken): version[0] == "5" returns False for "6.1.1"
# After (fixed): int(version.split('.')[0]) >= 5 returns True

# New version properties available:
orca = ade.methods.ORCA()
print(orca.is_v6)           # True for ORCA 6.x
print(orca.is_v5_or_later)  # True for ORCA 5.x and 6.x
print(orca.major_version)   # 6
```

## Dependencies

* [Python](https://www.python.org/) >= 3.9
* **Recommended High-Level Method:**
   * [ORCA](https://www.faccts.de/orca/) >= 6.0 (**Recommended** - free, full feature support)
* Alternative High-Level Methods:
   * [ORCA](https://www.faccts.de/orca/) 4.x-5.x (limited features)
   * [Gaussian09](https://gaussian.com/glossary/g09/)
   * [Gaussian16](https://gaussian.com/gaussian16/)
   * [NWChem](http://www.nwchem-sw.org/index.php/Main_Page) >= 6.5
   * [QChem](https://www.q-chem.com/) >= 5.4
   * [GPU4PySCF](https://github.com/pyscf/gpu4pyscf)
   * [TeraChem](https://www.petachem.com/products.html) >= 1.9 (GPU-accelerated, x86_64 only)
   * [CP2K](https://www.cp2k.org/) >= 2024.1 (GPU-accelerated molecular DFT)
* **Recommended Low-Level Method:**
   * [XTB](https://github.com/grimme-lab/xtb) >= 6.5 (**Recommended**)
* Alternative Low-Level Methods:
   * [g-xTB](https://github.com/grimme-lab/g-xtb) 2.0+ (Grimme's general-purpose semiempirical method, approximating ωB97M-V/def2-TZVPPD-quality results for Z = 1–103)
   * [MOPAC](http://openmopac.net/)

**Optional Dependencies:**
* [orca-pi](https://pypi.org/project/orca-pi/) - ORCA Python Interface for enhanced parsing
* [MLIP API Server](https://github.com/gomesgroup/mlip-api) - For MLIP-accelerated optimization

### GPU-Accelerated Wrappers

The gomesgroup fork includes wrappers for GPU-accelerated quantum chemistry codes:

| Code | Architecture | GPU Required | Use Case |
|------|--------------|--------------|----------|
| **TeraChem** | x86_64 only | NVIDIA (CUDA) | Fast GPU DFT, CASSCF, excited states |
| **CP2K** | x86_64 + ARM64 | NVIDIA (optional) | Molecular DFT with GPU acceleration |
| **GPU4PySCF** | x86_64 + ARM64 | NVIDIA (CC 7.0+) | Python-native GPU DFT |

**TeraChem Example:**
```python
import autode as ade
mol = ade.Molecule(smiles='c1ccccc1')
mol.optimise(method=ade.methods.TeraChem())
```

**CP2K Example (Molecular Systems):**
```python
import autode as ade
mol = ade.Molecule(smiles='O')
mol.single_point(method=ade.methods.CP2K())
mol.optimise(method=ade.methods.CP2K())
```

**Note:** CP2K wrapper is configured for molecular (non-periodic) calculations using WAVELET Poisson solver.

## Installation

To install **autodE** with [conda](https://anaconda.org/conda-forge/autode):
```
conda install autode -c conda-forge
```

To install the gomesgroup fork with ORCA 6.x features:
```
pip install git+https://github.com/gomesgroup/autode.git
```

See the [installation guide](https://duartegroup.github.io/autodE/install.html) for installing from source.

## Usage

Reaction profiles in **autodE** are generated by initialising _Reactant_ and _Product_ objects,
generating a _Reaction_ from those and invoking _calculate_reaction_profile()_.
For example, to calculate the profile for a 1,2 hydrogen shift in a propyl radical:

```python
import autode as ade

ade.Config.n_cores = 8
ade.Config.hcode = "orca"  # Recommended: ORCA 6.x
ade.Config.lcode = "xtb"

r = ade.Reactant(name='reactant', smiles='CC[C]([H])[H]')
p = ade.Product(name='product', smiles='C[C]([H])C')

reaction = ade.Reaction(r, p, name='1-2_shift')
reaction.calculate_reaction_profile()  # creates 1-2_shift/ and saves profile
```

### Using ORCA 6.x Features

**Modern DFT with D4 Dispersion:**
```python
import autode as ade
from autode.wrappers.keywords.orca6 import orca6_keywords, r2scan3c_keywords

ade.Config.hcode = "orca"

# Use r2SCAN-3c (fast, accurate composite method)
mol = ade.Molecule(smiles='c1ccccc1')
mol.optimise(method=ade.methods.ORCA(), keywords=r2scan3c_keywords)
```

**IRC Validation of Transition States:**
```python
from autode.wrappers import run_irc, validate_ts_with_irc

# Find TS
ts = reaction.find_ts()

# Validate with IRC (traces path both directions)
irc_result = validate_ts_with_irc(ts, method=ade.methods.ORCA())
if irc_result.connects_to_reactants and irc_result.connects_to_products:
    print("TS validated!")
```

**Hybrid MLIP-DFT NEB:**
```python
from autode.wrappers.neb import NEBMLIPHybrid

# Fast NEB with MLIP, refined with DFT at key points
neb = NEBMLIPHybrid(
    reactant=r, product=p,
    mlip_model='aimnet2',
    dft_method=ade.methods.ORCA()
)
ts_guess = neb.get_ts_guess()
```

See [examples/](https://github.com/duartegroup/autodE/tree/master/examples) for
more examples and [duartegroup.github.io/autodE/](https://duartegroup.github.io/autodE/) for
additional documentation.

## ORCA 6.x Feature Reference

| Feature | Module | Description |
|---------|--------|-------------|
| GOAT | `autode.wrappers.keywords.goat` | Metadynamics conformer search |
| SOLVATOR | `autode.wrappers.keywords.solvator` | Explicit solvation shell generation |
| DOCKER | `autode.wrappers.docker` | Molecular docking for complexes |
| IRC | `autode.wrappers.irc` | Intrinsic reaction coordinate |
| NEB-TS | `autode.wrappers.keywords.neb` | Enhanced nudged elastic band |
| MLIP | `autode.wrappers.mlip_external` | Machine learning potential integration |
| OPI | `autode.wrappers.opi_wrapper` | ORCA Python Interface |
| D4 | `autode.wrappers.keywords.dispersion` | D4 dispersion correction |
| r2SCAN-3c | `autode.wrappers.keywords.orca6` | Composite DFT method |
| Relativistic | `autode.wrappers.keywords.relativistic` | ZORA, X2C, DKH2 |

For detailed documentation on ORCA 6.x features, see [docs/orca6_features.md](docs/orca6_features.md).

## Multi-Architecture Support

The gomesgroup fork supports multiple architectures for MLIP integration:

| Architecture | Platform | MLIP Server |
|--------------|----------|-------------|
| x86_64 | Linux (Intel/AMD) | `localhost:5003` |
| aarch64 | Linux (ARM64/Grace Hopper) | `localhost:5003` |
| arm64 + Darwin | macOS (Apple Silicon) | `localhost:5003` |

## Development

There is a [slack workspace](https://autodeworkspace.slack.com) for development and discussion - please
[email](mailto:autodE-gh@outlook.com?subject=autodE%20slack) to be added. Pull requests are
very welcome but must pass all the unit tests prior to being merged. Please write code and tests!
See the [todo list](https://github.com/duartegroup/autodE/projects/1) for features on the horizon.
Bugs and feature requests should be raised on the [issue page](https://github.com/duartegroup/autodE/issues).

> **_NOTE:_**  We'd love more contributors to this project!

### Running Tests

```bash
# Run all tests
pytest tests/

# Run ORCA 6.x specific tests
pytest tests/test_orca_v6*.py -v

# Run with coverage
pytest tests/ --cov=autode --cov-report=html
```

## Citation

If **autodE** is used in a publication please consider citing the [paper](https://doi.org/10.1002/anie.202011941):

```
@article{autodE,
  doi = {10.1002/anie.202011941},
  url = {https://doi.org/10.1002/anie.202011941},
  year = {2021},
  publisher = {Wiley},
  volume = {60},
  number = {8},
  pages = {4266--4274},
  author = {Tom A. Young and Joseph J. Silcock and Alistair J. Sterling and Fernanda Duarte},
  title = {{autodE}: Automated Calculation of Reaction Energy Profiles -- Application to Organic and Organometallic Reactions},
  journal = {Angewandte Chemie International Edition}
}
```


## Contributors

- Tom Young ([@t-young31](https://github.com/t-young31))
- Joseph Silcock ([@josephsilcock](https://github.com/josephsilcock))
- Kjell Jorner ([@kjelljorner](https://github.com/kjelljorner))
- Thibault Lestang ([@tlestang](https://github.com/tlestang))
- Domen Pregeljc ([@dpregeljc](https://github.com/dpregeljc))
- Jonathon Vandezande ([@jevandezande](https://github.com/jevandezande))
- Shoubhik Maiti ([@shoubhikraj](https://github.com/shoubhikraj))
- Daniel Hollas ([@danielhollas](https://github.com/danielhollas))
- Nils Heunemann ([@nilsheunemann](https://github.com/NilsHeunemann))
- Sijie Fu ([@sijiefu](https://github.com/SijieFu))
- Javier Alfonso ([@javialra97](https://github.com/javialra97))
- Robert MacKnight ([@rmacknight99](https://github.com/rmacknight99))
- Jose Regio ([@jregio](https://github.com/jregio))
- Gabe Gomes ([@gabegomes](https://github.com/gabegomes)) - gomesgroup fork maintainer
