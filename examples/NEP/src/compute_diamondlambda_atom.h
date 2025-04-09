/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(diamondlambda/atom,ComputeDiamondLambdaAtom)

#else

#ifndef LMP_COMPUTE_DIAMONDLAMBDA_ATOM_H
#define LMP_COMPUTE_DIAMONDLAMBDA_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeDiamondLambdaAtom : public Compute {
 public:
  ComputeDiamondLambdaAtom(class LAMMPS *, int, char **);
  ~ComputeDiamondLambdaAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);

 private:
  int nmax,maxneigh,nnn,ndegree; 
  double cutsq,rsoft,cutbig;
  class NeighList *list;
  int **hydrogenBondNeigh_big;
  int **hydrogenBondNeigh;
  int *NearestNeighNumber;     
  double *distsqO,*distsqH;
  int *nearestO,*nearestH;
  int hydrogenId,oxygenId;
  double hydroDev;
  bool packQlm,packSolid,packNuclei;
  bool computeNucleiId;
  bool biggest;    
  int nConditions;
  int *value2index;
  bool *compareDirection;// true for greater than, false for less than
  double *threshold;
  double hardNeighbourDistance;
  int hardNeighbourCount;

  double **qlmarray;
  double *qnvector;
  double *nucleiID;
  double *isSolid;

  void add_qlm_complex(int m,double frr,double normx,double normy,double normz,double *u,double *v);
  void add_qn_complex(int, int, double, double*, double*);
  void calc_qn_trig(double, double, double&, double&);
  void select2(int, int, double *, int *);

  double polar_prefactor(int, int, double);
  double associated_legendre(int, int, double);
  double smearing(double) const;

  bool hydrogenBond(int i,int j,int ncountH) const;
  bool checkSolid(int ii) const;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute hexorder/atom requires a pair style be defined

Self-explantory.

W: More than one compute hexorder/atom

It is not efficient to use compute hexorder/atom more than once.

*/
