# NEP Example (FFS with NEP Potential)

This example shows how to use Forward Flux Sampling (FFS) with a Neural Network Potential (NEP) in LAMMPS.

---

## Compatibility Warning

This example is compatible only with LAMMPS version 2017 (or close).  
It is not compatible with newer versions such as LAMMPS 2024.

If run with LAMMPS 2024 or later, you may see an error like:

```
ERROR: Illegal pair_style command; nep doesn't require any parameter (../pair_NEP.cpp:117)
```

This is due to interface changes in `pair_style nep` in recent LAMMPS versions.

---

## NEP Source Code and Patch Notes

The NEP implementation is originally from:

https://github.com/brucefan1983/NEP_CPU/tree/main

To make it compatible with LAMMPS 2017, we applied the following patch to `pair_NEP.cpp:

We commented out the following conditional block:

```
#if LAMMPS_VERSION_NUMBER >= 20201130
  centroidstressflag = CENTROID_AVAIL;
#else
  centroidstressflag = 2;
#endif
```

- `centroidstressflag` is not available in LAMMPS 2017
- `CENTROID_AVAIL` was introduced in later versions


These changes are necessary for successful compilation on LAMMPS 2017.

---

## Required Source Code Files and Compilation

You need to copy the following files into your LAMMPS 2017 `src/` directory:

From `examples/NEP/src/`:

- `compute_biggest.cpp`, `compute_biggest.h`
- `compute_diamondlambda_atom.cpp`, `compute_diamondlambda_atom.h`
- `ffs.cpp`, `ffs.h`
- `lammps.cpp`, `main.cpp`
- `USER-NEP/` folder (entire folder with patched NEP)


Compile LAMMPS (2017 version) with:

```
make yes-USER-NEP
make mpi
```


## Examples in the Folder


- Diamond Nucleation @ 3700 K, 15 GPa
- Graphite Nucleation @ 3700 K, 15 GPa

Each example directory contains:
- `lammps.input` – LAMMPS input settings
- `ffs.input` – FFS parameters
- `trajectory.in.txt` – empty initial file
- `LiquidC.data` – liquid carbon configuration
- `C_2022_NEP3.txt` – NEP potential file
- `pool/` – stores generated `.xyz` files
- `run_lammps.sh` – job submission script



You can download the FFS output(e.g., FFS log files, configurations) for carbon crystallization here:

https://github.com/WanyuZhao/FFS-LAMMPS/releases/tag/carbon_NEP_FFS_output


Refer to `docs/manual.md` for more information about FFS setup and analysis.

## Reproduction Instruction

Once LAMMPS is compiled with `USER-NEP`, run the FFS simulations for either graphite or diamond:

```bash
cd examples/NEP/{graphite,diamond}
bash run_lammps.sh
```

---

## Notes

This example is preserved for reference and historical reproducibility.  
If you plan to use NEP potentials in modern LAMMPS (2024+), refer to:

https://github.com/brucefan1983/NEP_CPU

and follow their updated documentation.
