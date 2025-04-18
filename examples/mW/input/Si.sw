# Stillinger-Weber parameters for various elements and mixtures
# multiple entries can be added to this file, LAMMPS reads the ones it needs
# these entries are in LAMMPS "metal" units:
#   epsilon = eV; sigma = Angstroms
#   other quantities are unitless

# format of a single entry (one or more lines):
#   element 1, element 2, element 3, 
#   epsilon, sigma, a, lambda, gamma, costheta0, A, B, p, q, tol

# Here are the original parameters in metal units, for Silicon from:
#
# Stillinger and Weber,  Phys. Rev. B, v. 31, p. 5262, (1985)
#
#        epsilon      sigma   a     lambda gamma costheta0
#        A            B             p    q   tol

Si Si Si 0.268373607  2.3925  1.80  23.15  1.20  -0.333333333333
         7.049556277  0.6022245584  4.0  0.0 0.0

