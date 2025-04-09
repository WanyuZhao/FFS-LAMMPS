/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <math.h>
#include <string.h>
#include "compute_biggest.h"
#include "update.h"
#include "modify.h"
#include "atom.h"
#include "comm.h"
#include "group.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "library.h"


using namespace LAMMPS_NS;

#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PERATOM 8
#define INVOKED_LOCAL 16

/* ---------------------------------------------------------------------- */

ComputeBiggest::ComputeBiggest(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
//  if (narg != 4) {
//      error->all(FLERR,"Illegal compute biggest command");
//  }

  int iarg=3;
  rgRatio=-1;
  makegroup=0;

  if (strncmp(arg[iarg],"c_",2) == 0) {
      int iCompute=modify->find_compute(&arg[iarg][2]);
      if (iCompute<0) {
          error->all(FLERR,"Illegal compute nuclei/atom command");
      }
      compute_nuclei=modify->compute[iCompute];
      if (compute_nuclei->peratom_flag!=1) {
          error->all(FLERR,"Illegal compute nuclei/atom command");
      }
  }
  if (strcmp(arg[iarg+1],"rgRatio") == 0) {
	  rgRatio = force->numeric(FLERR,arg[iarg+2]);
	  if (rgRatio<=0) error->all(FLERR,"Illegal compute nuclei/atom command");
	  iarg=5;
  }  
  if (strcmp(arg[iarg+1],"groupBig") == 0) {
	  if (iarg+3 > narg) error->all(FLERR,"Illegal compute nuclei/atom command");
      int n = strlen(arg[iarg+2]) + 1;
      makegroup = 1;
      groupname = new char[n];
      strcpy(groupname,arg[iarg+2]);
  }
  else {
      error->all(FLERR,"Illegal compute nuclei/atom command");
  }

  vector_flag = 1;
  size_vector = 2;
  extscalar = 0;
  extvector = 0;

  vector = new double[2];
}

/* ---------------------------------------------------------------------- */

ComputeBiggest::~ComputeBiggest()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeBiggest::init()
{
}

/* ----------------------------------------------------------------------
   compute the radius-of-biggest tensor of group of atoms
   around center-of-mass cm
   must unwrap atoms to compute Rg tensor correctly
------------------------------------------------------------------------- */

void ComputeBiggest::compute_vector()
{
  invoked_vector = update->ntimestep;

  {
      if (!(compute_nuclei->invoked_flag & INVOKED_PERATOM)) {
          compute_nuclei->compute_peratom();
          compute_nuclei->invoked_flag |= INVOKED_PERATOM;
      }
  }

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  tagint *tag = atom->tag;
  int maxTagLocal=0;
  double **x = atom->x; 
  int maxTag;
  int *countLocal,*countGlobal,*xGlobal,*yGlobal,*zGlobal;
  int nprocs = comm->nprocs;
  double low[3], high[3], xy, yz, xz;
  int periodicity[3], boxChange;
  double lengthx,lengthy,lengthz;
  
  lammps_extract_box(lmp, low, high, &xy, &yz, &xz, periodicity, &boxChange);
  lengthx=high[0]-low[0];
  lengthy=high[1]-low[1];
  lengthz=high[2]-low[2];
      
  for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      if (tag[i]>maxTagLocal) {
          maxTagLocal=tag[i];
      }
  }
  MPI_Allreduce(&maxTagLocal,&maxTag,1,MPI_INT,MPI_MAX,world);
    
      
  memory->create(countLocal,1+maxTag,"biggest:countLocal");
  memory->create(countGlobal,1+maxTag,"biggest:countGlobal");

   
  for (int i=0;i<=maxTag;i++) {
      countLocal[i]=0;
  }

  double *vector_nuclei=compute_nuclei->vector_atom;

  
  for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      countLocal[(int)(vector_nuclei[i])]++;
  }
    
  MPI_Allreduce(countLocal,countGlobal,1+maxTag,MPI_INT,MPI_SUM,world);
  

  int maxCount=-1,iMax;
  int lastmax=maxTag+1;
  double rgMax;
  while (1){
      maxCount = -1;
	  for (int i = 1; i <= maxTag; i++) {
		  int current=countGlobal[i];
		  if (current>=lastmax) continue;
		  if (current>maxCount) {
			  maxCount=current;
			  iMax=i;
		  }  
	  }
          if (maxCount == lastmax){
		      maxCount=-1;
         	  break;
          }
	  if (rgRatio==-1) break;      
	  lastmax = maxCount;
	  int countClusterLocal = 0;
	  for (int i = 0; i < nlocal; i++) {
		if (!(mask[i] & groupbit)) continue;
		if ((int)(vector_nuclei[i])==iMax){
			countClusterLocal++;
		}
	  }
	  double *xLocal,*yLocal,*zLocal;
          memory->create(xLocal, countClusterLocal, "biggest:xLocal");
	  memory->create(yLocal, countClusterLocal, "biggest:yLocal");	
	  memory->create(zLocal, countClusterLocal, "biggest:zLocal");	 
	  int ncountClusterLocal = 0;
	  for (int i = 0; i < nlocal; i++) {
		if (!(mask[i] & groupbit)) continue;
		if ((int)(vector_nuclei[i])==iMax){
			xLocal[ncountClusterLocal]=x[i][0];
			yLocal[ncountClusterLocal]=x[i][1];
			zLocal[ncountClusterLocal]=x[i][2];
			ncountClusterLocal++;
		}
	  }
	  int *recvcounts, *displs;
	  double *xGlobal,*yGlobal,*zGlobal;
	  memory->create(recvcounts, nprocs, "biggest:recvcounts");
	  memory->create(displs, nprocs, "biggest:displs");
	  memory->create(xGlobal, maxCount, "biggest:xGlobal");
	  memory->create(yGlobal, maxCount, "biggest:yGlobal");	
	  memory->create(zGlobal, maxCount, "biggest:zGlobal");					
	  MPI_Allgather(&countClusterLocal, 1, MPI_INT, recvcounts, 1, MPI_INT, world);
	  displs[0] = 0;
	  for (int iproc = 1; iproc < nprocs; iproc++)
		 displs[iproc] = displs[iproc - 1] + recvcounts[iproc - 1];			
	  // allgatherv acquires list of populated IDs from all procs
	  MPI_Allgatherv(xLocal, countClusterLocal, MPI_DOUBLE, xGlobal, recvcounts, displs, MPI_DOUBLE, world);
	  MPI_Allgatherv(yLocal, countClusterLocal, MPI_DOUBLE, yGlobal, recvcounts, displs, MPI_DOUBLE, world);
	  MPI_Allgatherv(zLocal, countClusterLocal, MPI_DOUBLE, zGlobal, recvcounts, displs, MPI_DOUBLE, world);
	  double deltx,delty,deltz,deltxx=0.0,deltyy=0.0,deltzz=0.0;
	  for (int i = 0; i<maxCount; i++){
		  for (int j= i+1; j<maxCount; j++){
				deltx=abs(xGlobal[i]-xGlobal[j]);
				delty=abs(yGlobal[i]-yGlobal[j]);
				deltz=abs(zGlobal[i]-zGlobal[j]);
				deltx=deltx-floor(deltx/lengthx+0.5)*lengthx;
				delty=delty-floor(delty/lengthy+0.5)*lengthy;
				deltz=deltz-floor(deltz/lengthz+0.5)*lengthz;
				deltxx +=deltx*deltx;
				deltyy +=delty*delty;
				deltzz +=deltz*deltz;
			}
		}
	  memory->destroy(recvcounts);
	  memory->destroy(displs);
	  memory->destroy(xLocal);
	  memory->destroy(yLocal);
	  memory->destroy(zLocal);
	  memory->destroy(xGlobal);
	  memory->destroy(yGlobal);
	  memory->destroy(zGlobal);
	  double rg[6];
	  rg[0] = rg[1] = rg[2] = rg[3] = rg[4] = rg[5] = 0.0;  
	  rg[0] = sqrt(deltxx)/sqrt(deltyy);
	  rg[1] = sqrt(deltxx)/sqrt(deltzz);
	  rg[2] = sqrt(deltyy)/sqrt(deltzz);
	  rg[3] = sqrt(deltyy)/sqrt(deltxx);
	  rg[4] = sqrt(deltzz)/sqrt(deltxx);
	  rg[5] = sqrt(deltzz)/sqrt(deltyy);
	  rgMax=-1.0;
	  for (int j=0; j<=5; j++){
		  if (rg[j]>rgMax){
			  rgMax=rg[j];
		  }
	  }
	  if (rgMax<rgRatio) break;
  }
  
  if (makegroup == 1) {
	  int *flags;
	  memory->create(flags,nlocal,"biggest:flags");
      for (int i = 0; i < nlocal; i++){
			  flags[i]=0;
	}
      for (int i = 0; i < nlocal; i++){
          if (!(mask[i] & groupbit)) continue;
		  if ((int)(vector_nuclei[i])==iMax){
			  flags[i]=1;
		  }
	  }
      group->create(groupname,flags);
      memory->destroy(flags);
  }
  
  vector[0]=maxCount;
  vector[1]=rgMax;

  memory->destroy(countGlobal);
  memory->destroy(countLocal);
}
