# Custom Patch to LAMMPS for Forward Flux Sampling (FFS)

This repository contains a custom modification to the LAMMPS source code (`lammps.cpp` and `main.cpp`) from the official April 2024 release to enable a interface for Forward Flux Sampling (FFS).

I. Custom Modification 

1. `src/lammps.cpp`
The file was patched (line 465 to line 468) to invoke FFS by a new command-line option:


```
// ffs
else if (strcmp(arg[iarg], "-ffs") == 0) {
    iarg += 3;
}
```

2. `src/main.cpp`
The file was patched (line 71 to line 105) to enable launching FFS simulations from the main function.


```
#include "ffs.h"
...

if (!ffsRequested(argc,argv)) {

  // run LAMMPS normally
  
} else {
  ffs_main(argc,argv);  // run FFS
}
```


This change allows `main.cpp` to detect whether `-ffs` was called, and if so, skip normal LAMMPS execution and launch the FFS routine directly.

II. Compatibility
This is a custom extension. The `-ffs` flag is not supported in the official LAMMPS package.

