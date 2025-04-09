# Forward Flux Sampling with LAMMPS

This repository provides an integrated implementation of Forward Flux Sampling (FFS) within the LAMMPS. It includes both the core FFS code and a minimal custom patch to LAMMPS (April 2024 version).

---

## Installation

### Requirements
- LAMMPS version ≥ April 2024
- MPI-enabled system

### Compilation Steps
1. Download LAMMPS package 
2. Copy files in ffs-lammps/src to /lammps/src/ 
3. Compile LAMMPS normally. Follow the official instructions here: https://docs.lammps.org/Build.html


---

## Usage
See docs/manual.md for full details.

---

## Custom LAMMPS Patch

This repo includes minimal changes to:
- `lammps.cpp`: adds `-ffs` option
- `main.cpp`: routes execution to `ffs_main()`

See docs/custom_lammps_patch.md for full patch details.

---

## Examples

This repository contains two example systems located in the `examples/` directory:

- `examples/mW/` – corse-grained mW water system 
- `examples/NEP/` – Carbon system with NEP potential (diamond or graphite nucleation)

Each example includes:
- An `input/` folder containing the LAMMPS input files
- A corresponding `output.zip` archive containing the simulation output

---

## License

This code is under the terms of the GNU General Public License (GPL)(https://www.gnu.org/licenses/gpl-3.0.html), consistent with the licensing of the LAMMPS codebase, from which it is derived.

It includes modified versions of core LAMMPS files (`main.cpp`, `lammps.cpp`), and must therefore remain GPL-compliant.

---

## Image License
 
Workflow figure in this repository (e.g., `docs/workflow.png`) is original creation by the authors.

They are released under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

You are free to use these images in academic papers or presentations, provided appropriate attribution is given  
(e.g., "© 2025 Zhaowanyu, CC BY 4.0").

---


