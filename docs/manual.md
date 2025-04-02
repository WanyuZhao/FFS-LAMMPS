# FFS Manual

This document provides a complete overview of input files preparation, LAMMPS command usage, and output file postprocessing for Forward Flux Sampling (FFS) simulations.

---

##  Input Files Setup

Create a folder 'foldername'/ containing the following files to run FFS simulation:

```
'foldername'/
├── in.data              # Initial configuration of the liquid system
├── lammps.input         # Standard LAMMPS MD settings
├── ffs.input            # FFS parameters
├── Si.sw                # potential file (e.g., mW water)
├── trajectory.in.txt    # empty initial file; used by FFS to track progress
├── pool/                # empty folder; output configuration `.xyz` files will be automatically generated here
└── run_lammps.sh        # Job submission script to run FFS on HPC cluster
```

###  ffs.input

This file contains FFS parameters. Below is an example and description for each parameter (!!!!):

```
equilibrium 100000        # Equilibrium steps before FFS begins
check_every 20            # Frequency (in steps) to check lambda for every shooting
print_every 50            # Print every 50th lambda check result
water_group all           # Group name used to define the liquid system
temperature 220           # Target temperature of the simulation

config_each_lambda 120 120 120 120 120 120 120 120 120 120 120
                         # Number of successful configurations collected at interface

lambda 10 40 60 80 100 120 145 170 200 230 260 300
                         # Lambda value at each interface
```


### lammps.input

Below is a description of the usage, syntax, arguments, and examples for the compute command:

####  compute diamondlambda/atom

This compute identifies the crystalline-like molecules by employing the local order parameter, and then determines clusters consisted with crystalline-like molecules.  
This compute returns tag ID of each identified cluster which will be used in compute biggest command.

Syntax:

compute ID group_ID diamondlambda/atom degree nnn cutoff cutoff_big orderParameterOnly nucleiBiggest self argus


Arguments:
- `ID`: User-defined name assigned to identify this compute command
- `group_ID`: Specifies the ID of the group of atoms
- `diamondlambda/atom`: Specifies the style name of this compute command
- `degree`: Degree of local order parameter to identify the crystalline cluster in liquid (e.g., 6 for ice) 
- `nnn`: Required number of nearest neighbors crystalline-like cluster must have
- `cutoff`: Distance cutoff for identifying nearest neighbors
- `cutoff_big`: Cutoff used for cluster connectivity
- `orderParameterOnly`: Return only local order parameters
- `nucleiBiggest`: Counts first-neighbor molecules as part of the crystalline cluster even if they do not individually satisfy the local order parameter threshold
- `self`: Defines the direction of the local order parameter threshold. Must be followed by greaterThan or lessThan to select molecules with values above or below the threshold, respectively

Examples:

compute iceId water diamondlambda/atom degree 6 nnn 4 cutoff 3.2 cutoff_big 3.2 nucleiBiggest self greaterThan 0.5
compute graphiteId graphite diamondlambda/atom degree 3 nnn 3 cutoff 1.8 orderParameterOnly




####  compute biggest

This compute identifies the biggest crystalline-like cluster and return its atoms number.

Syntax:

compute lambda group_ID biggest c_ID groupBig cluster_ID


Arguments:
- `lambda`: ID for this compute command; must always be set to lambda
- `group_ID`: Specifies the ID of the group of atoms
- `biggest`: Specifies the style name of this compute command
- `c_ID`: The compute ID from `diamondlambda/atom`
- `groupBig`:Defines a group for the atoms in the biggest crystalline-like cluster; followed by `cluster_ID`, which assigns a name to this group

Example:
compute lambda water biggest c_iceId groupBig biggestcluster

---

## Running the Simulation


mpirun -np 512 ./lmp_mpi -in lammps.input -screen none -ffs 128 ffs.input

- `-ffs` enables FFS mode
- `128`: Number of MD shootings per interface

---


## Output Files

### trajectory.out.txt

Each line corresponds to a successful trial that reaches the next interface.

Example line:

```
49 48 (xyz.0__40_0) >==2016267247 400==> 63 (xyz.1__49_0)
```

This line corresponding to a successful trail reaches lambda_1 from lambda_0

Fields:
1. `49`: Index of the shooting
2. `48`: Size of the biggest crystalline cluster in the parent config
3. `xyz.0__40_0`: Name of the parent config (interface = 0, shooting index = 40, 0 = sequence number of the successful trial resulting from shooting 40)
4. `2016267247`: Random seed used for reproducibility
5. `400`: Timesteps to reach next interface
6. `63`: Size of biggest crystalline cluster in the child config
7. `xyz.1__49_0`: Name of child configuration (interface = 1, shooting index = 49, 0 = sequence number of the successful trial resulting from shooting 49)


### slurm-xxxx.out

This file includes:
- Status printouts during `check_every` and `print_every`
- Lines indicating failed and successful attempts
- MPI communication signals for internal bookkeeping

Example status printout line:

```
[date=1723577788] [universe=38] [steps=100] : 8 ... 40
```
Information output for every lambda check. Associated with `check_every` and `print_every` command in the ffs.input.

Fields:
2. `8`: size of largest crystalline cluster of shooting(universe) 38
3. `40`: Target size of the next interface



Example failure line:
```
48 (xyz.0__10_2) >==1490733239 40==> 11 (__________)
```
This indicates the trail returned from lambda_0 to lambda_A.

Fields:
1. `48`: Size of the biggest crystalline cluster in the parent config
2. `xyz.0__10_2`: Name of the parent config (interface = 0, shooting index = 10, 2 = sequence number of the successful trial resulting from shooting 10)
4. `1490733239`: Random seed used for reproducibility
5. `40`: Timesteps to reach next interface

Example successful line is the same as example line in the trajectory.out.txt



MPI Tags:
```
TAG_FILEWRITER_LINE
TAG_COUNTDOWN_DONE
```
These confirm communication between universes and the world leader.


---

## Postprocessing

Use the included Python script in `examples/` to compute the nucleation rate:
```
python NucleationRate.py --datafile in.data --logfile slurm-xxxx.out
```
