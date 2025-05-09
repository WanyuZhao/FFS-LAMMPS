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

- **examples/mW/**: coarse-grained mW water system  
  - `input/` — LAMMPS input files  
  - `output.zip` — FFS simulation outputs
  
- **examples/NEP/**: carbon system with NEP potential (diamond or graphite nucleation)  
  - LAMMPS input files are provided in the `examples/NEP/graphite/` and `examples/NEP/diamond/` folders  
  - FFS outputs can be downloaded from the [carbon_NEP_FFS_output release](https://github.com/WanyuZhao/FFS-LAMMPS/releases/tag/carbon_NEP_FFS_output)    


---

## License


This package (including all patches to LAMMPS core files such as `main.cpp` and `lammps.cpp`, documentation, and original workflow figures) is licensed under the GNU General Public License version 2.0 (GPL-2.0). See the LICENSE file for full details.

This is consistent with the licensing of the LAMMPS codebase from which it is derived.


---



