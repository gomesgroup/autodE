# autodE ORCA 6.x Features Guide

This document covers the new ORCA 6.x features integrated into autodE, including installation instructions, usage examples, migration guidance, and API reference.

**Document Version:** January 2026
**autodE Version:** 1.4.4+ (with ORCA 6.x extensions)
**Supported ORCA Versions:** 5.x, 6.0, 6.1

---

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [New Features](#new-features)
   - [GOAT Workflow](#goat-workflow)
   - [SOLVATOR (Explicit Solvation)](#solvator-explicit-solvation)
   - [DOCKER (Molecular Docking)](#docker-molecular-docking)
   - [IRC (Intrinsic Reaction Coordinate)](#irc-intrinsic-reaction-coordinate)
   - [NEB-TS (Updated Syntax)](#neb-ts-updated-syntax)
   - [MLIP Acceleration via ExtOpt](#mlip-acceleration-via-extopt)
4. [Migration Guide](#migration-guide)
5. [API Reference](#api-reference)
6. [Troubleshooting](#troubleshooting)

---

## Overview

ORCA 6.x introduces several major features that enhance computational chemistry workflows:

| Feature | Description | ORCA Version |
|---------|-------------|--------------|
| **GOAT** | Global Optimization And Transition state search | 6.0+ |
| **SOLVATOR** | Automated explicit solvation shell generation | 6.0+ |
| **DOCKER** | QM-level molecular docking with GFNn-xTB | 6.0+ |
| **IRC** | Intrinsic Reaction Coordinate with new algorithms | 6.0+ |
| **NEB-TS** | Enhanced NEB with NEB_TS_XYZ input format | 6.0+ |
| **D4 Dispersion** | New default dispersion correction | 6.0+ |
| **ExtOpt** | External optimizer interface for MLIP | 6.0+ |

### Version Compatibility

| ORCA Version | D3BJ Default | D4 Default | GOAT | SOLVATOR | DOCKER | NEB_TS_XYZ |
|--------------|--------------|------------|------|----------|--------|------------|
| 5.0.x | Yes | No | No | No | No | No |
| 6.0.x | No | Yes | Yes | Yes | Yes | Yes |
| 6.1.x | No | Yes | Yes | Yes | Yes | Yes |

autodE automatically detects your ORCA version and adjusts keyword defaults accordingly.

---

## Installation

### Prerequisites

- ORCA 6.x installed and available in PATH
- Python 3.9+
- autodE 1.4.4+

### Basic Installation

```bash
# Install autodE with all dependencies
pip install autode

# Or on the GPG cluster:
source /mnt/beegfs/software/autode/setup-autode.sh
```

### OPI (ORCA Python Interface) Installation

OPI is the official Python interface for ORCA 6.x, providing type-safe input generation and structured output parsing.

```bash
# Install OPI (optional but recommended for ORCA 6.x)
pip install orca-pi
```

**Benefits of OPI:**
- Type-safe input file generation
- Structured JSON output parsing
- Access to ORCA 6.x-specific features
- Better error handling and validation

### orca-external-tools Setup

For MLIP acceleration via ExtOpt, you need the client script:

```bash
# On the GPG cluster, this is pre-configured:
source /mnt/beegfs/software/autode/setup-autode.sh

# The client script is at:
# /mnt/beegfs/software/autode/bin/mlip_orca_client.sh
```

### Multi-Architecture Notes

autodE with ORCA 6.x works on:

| Architecture | Nodes | Status | Notes |
|--------------|-------|--------|-------|
| **x86_64** | gpg-boltzmann, EPYC nodes | Full support | All features including MLIP |
| **ARM64** (aarch64) | Grace Hopper (pearl, kahneman, gh-3, gh-4) | Full support | MLIP via remote server |
| **Apple Silicon** | macOS M1/M2/M3 | Partial | MLIP requires network access to server |

**Important:** The MLIP server runs on gpg-boltzmann:5003. ARM64 and Apple Silicon clients connect remotely.

---

## New Features

### GOAT Workflow

GOAT (Global Optimization And Transition state search) uses metadynamics-based exploration to find:
- Global minima (conformer search)
- Transition states
- Reaction intermediates

#### GOAT Calculation Types

| Type | Keyword | Description |
|------|---------|-------------|
| `OPT` | `GOAT_OPT` | Find global minimum via conformer search |
| `TS` | `GOAT_TS` | Find transition states via NEB-TS or rmsPATH |
| `REACT` | `GOAT_REACT` | Explore reaction space for multi-molecule systems |

#### Basic Usage

```python
from autode.wrappers.keywords import GOATKeywords, goat_conformer_search

# Create GOAT keywords for conformer search
goat_kws = GOATKeywords(
    goat_type="OPT",
    method="r2SCAN-3c",
    ewin=20.0,      # Energy window in kcal/mol
    max_opt=100,    # Maximum optimizations
    md_length=250,  # Metadynamics length in fs
)

# Generate ORCA input
print(goat_kws.to_orca_input())
```

**Output:**
```
!GOAT_OPT r2SCAN-3c
%goat
  Ewin 20.0
  MaxOpt 100
  MD_Length 250
end
```

#### Conformer Search with GOAT

```python
from autode import Molecule
from autode.wrappers.keywords import goat_conformer_search

# Load molecule
mol = Molecule(smiles="CCCCCC")  # Hexane

# Run GOAT conformer search
conformers = goat_conformer_search(
    molecule=mol,
    method="r2SCAN-3c",
    ewin=15.0,
    max_opt=50,
    n_cores=8,
)

print(f"Found {len(conformers)} unique conformers")
```

#### TS Search with GOAT

```python
from autode import Molecule
from autode.wrappers.keywords import goat_ts_search

# Load reactant
reactant = Molecule("reactant.xyz")

# Find transition states
result = goat_ts_search(
    reactant=reactant,
    method="r2SCAN-3c",
    ts_search_method="NEB-TS",  # or "rmsPATH"
    ewin=20.0,
    n_cores=8,
)

if result['ts_candidates']:
    print(f"Found {len(result['ts_candidates'])} TS candidates")
    best_ts = result['ts_candidates'][0]
    print(f"Best TS energy: {best_ts['energy']:.6f} Eh")
    print(f"Imaginary frequency: {best_ts['imaginary_freq']:.1f} cm^-1")
```

#### Parsing GOAT Output

```python
from autode.wrappers.keywords import parse_goat_output, parse_goat_allxyz

# Parse ORCA output file
results = parse_goat_output("goat_calculation.out")
print(f"Calculation type: {results['calculation_type']}")
print(f"Found {results['n_conformers']} conformers")
print(f"Terminated normally: {results['terminated_normally']}")

# Parse .allxyz trajectory file
conformers = parse_goat_allxyz("goat_calculation.GOAT.allxyz")
for conf in conformers:
    print(f"Conformer {conf['index']}: E = {conf['energy']:.6f} Eh")
```

---

### SOLVATOR (Explicit Solvation)

SOLVATOR automatically generates explicit solvation shells around a solute using force field-based placement.

#### Supported Solvents

| Solvent | Aliases |
|---------|---------|
| water | h2o |
| methanol | - |
| ethanol | - |
| acetonitrile | mecn, acn |
| dmso | dimethylsulfoxide |
| thf | tetrahydrofuran |
| acetone | - |
| dichloromethane | dcm, ch2cl2, methylenechloride |
| chloroform | chcl3 |
| toluene | - |
| benzene | - |
| hexane | n-hexane |
| diethylether | ether, et2o |
| dioxane | - |
| pyridine | - |
| dmf | dimethylformamide |

#### Basic Usage

```python
from autode.wrappers.solvator import SolvatorKeywords, solvate_molecule

# Create SOLVATOR keywords
solv_kws = SolvatorKeywords(
    solvent="water",
    n_shells=2,           # Number of solvation shells
    force_field="GFN-FF"  # Force field for placement
)

# Generate ORCA input
print(solv_kws.to_orca_keyword())  # "SOLVATOR"
print(solv_kws.to_orca_block())
```

**Output:**
```
%solvator
  Solvent "water"
  SolvShell 2
  ForceField "GFN-FF"
end
```

#### Complete SOLVATOR Input

```python
from autode.wrappers.solvator import generate_solvator_input, SolvatorKeywords

# Define molecule coordinates
coordinates = [
    ("O", 0.0, 0.0, 0.1173),
    ("H", 0.0, 0.7572, -0.4692),
    ("H", 0.0, -0.7572, -0.4692),
]

# Generate complete ORCA input
input_text = generate_solvator_input(
    keywords=SolvatorKeywords(solvent="water", n_shells=2),
    charge=0,
    multiplicity=1,
    coordinates=coordinates,
    method="B3LYP",
    basis="def2-SVP",
)

print(input_text)
```

#### Advanced Configuration

```python
from autode.wrappers.solvator import SolvatorKeywords, SolvatorConfig

# Advanced configuration with counterions
config = SolvatorConfig(
    counterion="Cl-",
    n_counterions=2,
    seed=42,
    max_iterations=2000,
)

solv_kws = SolvatorKeywords(solvent="water", n_shells=3)
print(solv_kws.to_orca_block(config=config))
```

**Output:**
```
%solvator
  Solvent "water"
  SolvShell 3
  ForceField "GFN-FF"
  Counterion "Cl-"
  NCounterions 2
  Seed 42
  MaxIterations 2000
end
```

---

### DOCKER (Molecular Docking)

DOCKER performs QM-level molecular docking to generate binding poses for host-guest systems.

#### Basic Usage

```python
from autode.wrappers.docker import dock_molecules, DockerPose
from autode import Molecule

# Load host and guest molecules
host = Molecule("receptor.xyz")
guest = Molecule("ligand.xyz")

# Run DOCKER docking
poses = dock_molecules(
    host=host,
    guest=guest,
    n_poses=10,
    ewin=10.0,        # Energy window in kcal/mol
    gfn_method=2,     # GFN2-xTB
)

# Get best pose
if poses:
    best = poses[0]  # Sorted by binding energy
    print(f"Best binding energy: {best.binding_energy:.2f} kcal/mol")
```

#### DockerKeywords Class

```python
from autode.wrappers.keywords import DockerKeywords

# Create DOCKER keywords
docker_kws = DockerKeywords(
    guest_file="guest.xyz",
    gfn_method=2,
    poses=10,
    ewin=10.0,
    anchor_guest=0,   # Anchor atom index in guest
    anchor_host=5,    # Anchor atom index in host
)

print(docker_kws.to_orca_block())
```

**Output:**
```
%docker
  GuestFile "guest.xyz"
  GFN 2
  Poses 10
  Ewin 10.0
  AnchorGuest 0
  AnchorHost 5
end
```

#### Parsing DOCKER Output

```python
from autode.wrappers.docker import parse_docker_allxyz, filter_poses_by_ewin

# Parse .allxyz output file
poses = parse_docker_allxyz("docking.DOCKER.allxyz")

# Filter by energy window
filtered = filter_poses_by_ewin(poses, ewin=5.0)
print(f"Poses within 5 kcal/mol of best: {len(filtered)}")
```

#### Pre-Reactive Complex Generation

```python
from autode.wrappers.docker import generate_prereactive_complex
from autode import Molecule

# Generate pre-reactive complex for bimolecular reaction
reactant1 = Molecule("A.xyz")
reactant2 = Molecule("B.xyz")

complex = generate_prereactive_complex(
    reactant1=reactant1,
    reactant2=reactant2,
    n_poses=10,
    ewin=10.0,
)

print(f"Complex has {complex.n_atoms} atoms")
```

---

### IRC (Intrinsic Reaction Coordinate)

IRC traces the minimum energy path from a transition state to both reactants and products.

#### IRC Algorithms

| Algorithm | Description | Best For |
|-----------|-------------|----------|
| `LQA` | Local Quadratic Approximation | Most reactions (default) |
| `EULER` | Simple Euler stepping | Fast exploration |
| `MN` | Morokuma-Newton | High accuracy |

#### Basic Usage

```python
from autode.wrappers.keywords import IRCKeywords

# Create IRC keywords
irc_kws = IRCKeywords(
    direction="BOTH",    # FORWARD, BACKWARD, or BOTH
    max_iter=100,
    step_size=0.1,       # Bohr*sqrt(amu)
    algorithm="LQA",
)

# Add method keywords
irc_kws.append("r2SCAN-3c")

# Generate ORCA block
print(irc_kws.to_orca_keywords())
print(irc_kws.to_orca_block())
```

**Output:**
```
IRC r2SCAN-3c
%irc
  Direction BOTH
  MaxIter 100
  StepSize 0.1
  Algorithm LQA
end
```

#### TS Validation with IRC

```python
from autode.wrappers.irc import run_irc, validate_ts_with_irc, IRCPath
from autode import Molecule

# Load transition state
ts = Molecule("ts.xyz")

# Run IRC to validate TS
result = validate_ts_with_irc(
    ts_mol=ts,
    expected_reactant=Molecule("reactant.xyz"),
    expected_product=Molecule("product.xyz"),
    rmsd_threshold=0.5,  # Angstrom
    direction="BOTH",
    max_iter=100,
)

if result.is_valid:
    print("TS connects expected reactant and product!")
    print(f"Forward barrier: {result.irc_path.forward_barrier * 627.5:.1f} kcal/mol")
    print(f"Reverse barrier: {result.irc_path.reverse_barrier * 627.5:.1f} kcal/mol")
```

#### Parsing IRC Output

```python
from autode.wrappers.irc import IRCPath

# Parse IRC results from ORCA output
irc_path = IRCPath.from_orca_output("irc_calculation.out")

print(f"Forward steps: {irc_path.forward_steps}")
print(f"Reverse steps: {irc_path.reverse_steps}")
print(f"Converged: {irc_path.converged}")

# Get energies along the path
for point in irc_path:
    print(f"Step {point.step}: E = {point.energy:.6f} Eh ({point.direction})")
```

---

### NEB-TS (Updated Syntax)

ORCA 6.x introduces the `NEB_TS_XYZ` input format for direct XYZ file input.

#### NEBKeywords Class

```python
from autode.wrappers.keywords import NEBKeywords, NEBTSKeywords

# Create NEB keywords
neb_kws = NEBKeywords(
    n_images=12,
    spring=0.01,           # Eh/Bohr^2
    ts_search_algo="EF",   # EF or P-RFO
    interpolation="IDPP",  # Linear or IDPP
    free_end=False,
    preopt=True,
    max_iter=500,
)

print(neb_kws.to_block())
```

**Output:**
```
%neb
  NImages 12
  Spring 0.0100
  TS_Search_Algo EF
  Interpolation IDPP
  Free_End false
  PreOpt true
  MaxIter 500
end
```

#### NEB_TS_XYZ Input Format (ORCA 6.x)

```python
# Generate NEB_TS_XYZ coordinate block
print(neb_kws.to_neb_ts_xyz_block(
    reactant_filename="reactant.xyz",
    product_filename="product.xyz",
    charge=0,
    mult=1,
))
```

**Output:**
```
*NEB_TS_XYZ 0 1
reactant.xyz
product.xyz
```

#### Complete NEB-TS Input

```python
from autode.wrappers.keywords import NEBTSKeywords, NEBKeywords

# Create complete NEB-TS specification
nebts = NEBTSKeywords(
    neb_keywords=NEBKeywords(n_images=16, interpolation="IDPP"),
    method="r2SCAN-3c",
    additional_keywords=["TightSCF"],
)

print(nebts.generate_method_line())
# Output: !NEB-TS r2SCAN-3c TightSCF
```

#### Pre-defined NEB-TS Keyword Sets

```python
from autode.wrappers.keywords.orca6 import get_neb_ts_keywords

# Get customized NEB-TS keywords
kw_set = get_neb_ts_keywords(
    n_images=16,
    method="PBE0 def2-SVP D4",
    ts_search_algo="EF",
    interpolation="IDPP",
)
```

---

### MLIP Acceleration via ExtOpt

Use Machine Learning Interatomic Potentials (MLIP) to accelerate geometry optimization.

#### MLIPConfig Class

```python
from autode.wrappers.mlip_external import MLIPConfig

# Create MLIP configuration
config = MLIPConfig(
    server_url="http://gpg-boltzmann:5003",  # Auto-detected if None
    model="aimnet2",     # or "uma" for transition metals
    timeout=30,
    fallback_enabled=True,
)

# Check server availability
if config.is_server_available():
    print("MLIP server is available")
    models = config.get_available_models()
    print(f"Available models: {models}")
```

#### Supported MLIP Models

| Model | Alias | Description | Elements |
|-------|-------|-------------|----------|
| `aimnet+base` | `aimnet2` | AIMNet2 base model | H, C, N, O, F, S, Cl |
| `aimnet+pd` | - | AIMNet2 with periodic disturbance | Same |
| `aimnet+spin` | - | AIMNet2 with spin support | Same |
| `omol+uma_sm` | `uma` | Universal MLIP (includes transition metals) | Wide coverage |
| `omol+esen_sm_conserving` | - | Energy-conserving model | Wide coverage |

#### MLIPOptKeywords Class

```python
from autode.wrappers.mlip_external import MLIPOptKeywords, MLIPConfig

# Create MLIP-accelerated optimization keywords
config = MLIPConfig(model="aimnet2")
mlip_kws = MLIPOptKeywords(
    keyword_list=["def2-SVP"],  # Basis for SCF guess
    mlip_config=config,
)

# Generate ORCA %extopt block
print(mlip_kws.get_extopt_block())
```

**Output:**
```
%extopt
  CMD "/mnt/beegfs/software/autode/bin/mlip_orca_client.sh http://gpg-boltzmann:5003 aimnet+base"
end
```

#### Multi-Architecture Support

```python
from autode.wrappers.mlip_external import get_platform_info, get_default_server_url

# Check current platform
info = get_platform_info()
print(f"Architecture: {info['machine']}")
print(f"Is ARM64: {info['is_arm64']}")
print(f"Is Apple Silicon: {info['is_apple_silicon']}")

# Get appropriate server URL
url = get_default_server_url()
print(f"Default MLIP server: {url}")
```

---

## Migration Guide

### Deprecated Keywords and Replacements

| ORCA 5.x | ORCA 6.x | Notes |
|----------|----------|-------|
| `D3BJ` | `D4` | D4 is now default; D3BJ still works |
| `NEB ... *xyzfile` | `*NEB_TS_XYZ` | New input format for endpoints |
| `%neb MaxIter` | `%neb MaxIter` | Same (no change) |
| - | `GOAT_OPT` | New conformer search method |
| - | `SOLVATOR` | New explicit solvation |
| - | `DOCKER` | New molecular docking |

### D4 as New Default Dispersion

ORCA 6.x uses D4 dispersion by default. autodE automatically handles this:

```python
from autode.wrappers.keywords.orca6 import orca6_keywords, orca5_keywords

# ORCA 6.x keyword set (uses D4)
kws_v6 = orca6_keywords.opt  # Includes D4

# ORCA 5.x compatible (uses D3BJ)
kws_v5 = orca5_keywords.opt  # Includes D3BJ

# Or explicitly select version
from autode.wrappers.keywords.orca6 import get_orca_keywords
kws = get_orca_keywords(version="6")  # D4 default
kws = get_orca_keywords(version="5")  # D3BJ default
```

### Updating Existing Scripts

**Before (ORCA 5.x):**
```python
from autode import Molecule
from autode.methods import ORCA

mol = Molecule(smiles="CCCO")
mol.optimise(method=ORCA())  # Uses D3BJ
```

**After (ORCA 6.x):**
```python
from autode import Molecule
from autode.methods import ORCA

mol = Molecule(smiles="CCCO")
mol.optimise(method=ORCA())  # Automatically uses D4
```

No code changes needed - autodE detects ORCA version automatically!

---

## API Reference

### Classes

#### GOATKeywords

```python
class GOATKeywords(Keywords):
    """Keywords for ORCA GOAT calculations."""

    def __init__(
        self,
        goat_type: str,           # "OPT", "TS", or "REACT"
        method: str,              # e.g., "r2SCAN-3c"
        ewin: float = 20.0,       # Energy window (kcal/mol)
        max_opt: int = 100,       # Maximum optimizations
        md_length: int = 250,     # Metadynamics length (fs)
        ts_search_method: str = "NEB-TS",  # For GOAT_TS
        n_conformers: int = None, # Target conformer count
        keyword_list: list = None,
    ): ...

    def to_orca_input(self) -> str: ...
```

#### SolvatorKeywords

```python
class SolvatorKeywords:
    """Keywords for ORCA SOLVATOR explicit solvation."""

    def __init__(
        self,
        solvent: str = "water",
        n_shells: int = 2,
        force_field: str = "GFN-FF",
    ): ...

    def to_orca_keyword(self) -> str: ...
    def to_orca_block(self, config: SolvatorConfig = None) -> str: ...
```

#### DockerKeywords

```python
class DockerKeywords:
    """Keywords for ORCA DOCKER molecular docking."""

    def __init__(
        self,
        guest_file: str,
        gfn_method: int = 2,
        poses: int = 10,
        ewin: float = 10.0,
        anchor_guest: int = None,
        anchor_host: int = None,
    ): ...

    def to_orca_block(self) -> str: ...
```

#### IRCKeywords

```python
class IRCKeywords(Keywords):
    """Keywords for IRC calculations."""

    def __init__(
        self,
        keyword_list: list = None,
        direction: str = "BOTH",    # FORWARD, BACKWARD, BOTH
        max_iter: int = 100,
        step_size: float = 0.1,     # Bohr*sqrt(amu)
        algorithm: str = "LQA",     # LQA, EULER, MN
    ): ...

    def to_orca_block(self) -> str: ...
    def to_orca_keywords(self) -> str: ...
```

#### NEBKeywords

```python
class NEBKeywords:
    """Keywords for ORCA NEB-TS calculations."""

    def __init__(
        self,
        n_images: int = 8,
        spring: float = 0.01,
        ts_search_algo: str = "EF",      # EF or P-RFO
        interpolation: str = "IDPP",     # Linear or IDPP
        free_end: bool = False,
        preopt: bool = True,
        max_iter: int = 500,
    ): ...

    def to_block(self) -> str: ...
    def to_neb_ts_xyz_block(self, reactant_filename, product_filename, charge, mult) -> str: ...
```

#### MLIPConfig

```python
class MLIPConfig:
    """Configuration for MLIP server connection."""

    def __init__(
        self,
        server_url: str = None,    # Auto-detected if None
        model: str = "aimnet2",
        client_script: str = None,
        timeout: int = 30,
        fallback_enabled: bool = True,
        charge: int = 0,
        mult: int = 1,
    ): ...

    def is_server_available(self) -> bool: ...
    def get_available_models(self) -> list: ...
```

### Functions

| Function | Description |
|----------|-------------|
| `parse_goat_output(file)` | Parse GOAT ORCA output file |
| `parse_goat_allxyz(file)` | Parse GOAT .allxyz trajectory |
| `goat_workflow(mol, ...)` | Run complete GOAT workflow |
| `goat_conformer_search(mol, ...)` | GOAT-based conformer search |
| `goat_ts_search(reactant, ...)` | GOAT-based TS search |
| `generate_solvator_input(...)` | Generate SOLVATOR input file |
| `solvate_molecule(mol, ...)` | Solvate a molecule |
| `dock_molecules(host, guest, ...)` | Run DOCKER docking |
| `parse_docker_allxyz(file)` | Parse DOCKER output |
| `run_irc(ts_mol, ...)` | Run IRC calculation |
| `validate_ts_with_irc(ts, ...)` | Validate TS via IRC |

---

## Troubleshooting

### Common Issues

#### "GOAT_OPT not recognized"

**Problem:** ORCA reports unknown keyword `GOAT_OPT`

**Solution:** Ensure you have ORCA 6.0 or later installed:
```bash
orca --version  # Should show 6.x
```

#### "MLIP server not available"

**Problem:** `MLIPServerUnavailableError` when using MLIP acceleration

**Solution:**
1. Check if running on the cluster with network access
2. Verify server is running: `curl http://gpg-boltzmann:5003/models`
3. Set fallback enabled: `MLIPConfig(fallback_enabled=True)`

#### "D4 dispersion not found"

**Problem:** Error with D4 dispersion correction

**Solution:** Update to ORCA 6.x or use D3BJ:
```python
from autode.wrappers.keywords import d3bj
# Use d3bj instead of d4 for ORCA 5.x
```

#### "NEB_TS_XYZ syntax error"

**Problem:** ORCA doesn't recognize `*NEB_TS_XYZ` format

**Solution:** This is an ORCA 6.x feature. For ORCA 5.x, use the standard `*xyzfile` format:
```python
# ORCA 5.x compatible
# *xyzfile charge mult reactant.xyz
# product: product.xyz
```

### Getting Help

- autodE documentation: https://autodE.readthedocs.io/
- ORCA manual: https://orca-manual.mpi-muelheim.mpg.de/
- OPI documentation: https://github.com/faccts/opi
- GPG cluster issues: Contact cluster administrator

---

## References

1. ORCA 6.0 Manual - Max-Planck-Institut fur Kohlenforschung
2. autodE: A Python package for automated reaction path computation - T. A. Habershon et al.
3. OPI (ORCA Python Interface) - https://github.com/faccts/opi
4. D4 Dispersion Correction - E. Caldeweyher et al., J. Chem. Phys. 150, 154122 (2019)
5. r2SCAN-3c Composite Method - S. Grimme et al., J. Chem. Phys. 154, 064103 (2021)
