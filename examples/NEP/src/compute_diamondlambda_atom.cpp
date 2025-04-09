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

/* ----------------------------------------------------------------------
   Contributing author:  Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include <complex>
#include <string.h>
#include <stdlib.h>
#include "compute_diamondlambda_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeDiamondLambdaAtom::ComputeDiamondLambdaAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3 ) error->all(FLERR,"Illegal compute diamondlambda/atom command");

  ndegree = 6;
  nnn = 6;
  cutsq = 0.0;
  cutbig = 0.0;
  rsoft=0;
  hydroDev=-1;
  oxygenId=-1;
  hydrogenId=-1;
  computeNucleiId=true;
  biggest=false;
  
  nConditions=1;//This should be able to be extended in the future
  compareDirection=new bool[nConditions];
  threshold=new double[nConditions];
  hardNeighbourDistance=0;
  hardNeighbourCount=0;
  int iCondition=0;


  // process optional args

  int iarg = 3;
  while (iarg < narg) {
      if (strcmp(arg[iarg],"degree") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal compute diamondlambda/atom command");
          ndegree = force->numeric(FLERR,arg[iarg+1]);
          if (ndegree < 0)
              error->all(FLERR,"Illegal compute diamondlambda/atom command");
          iarg += 2;
      } else if (strcmp(arg[iarg],"rsoft") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal compute diamondlambda/atom command");
          rsoft = force->numeric(FLERR,arg[iarg+1]);
          if (rsoft <= 0)
              error->all(FLERR,"Illegal compute diamondlambda/atom command");
          iarg += 2;
      } else if (strcmp(arg[iarg],"nnn") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal compute diamondlambda/atom command");
          nnn = force->numeric(FLERR,arg[iarg+1]);
          if (nnn <= 0)
              error->all(FLERR,"Illegal compute diamondlambda/atom command");
          iarg += 2;
      } else if (strcmp(arg[iarg],"cutoff") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal compute diamondlambda/atom command");
          double cutoff = force->numeric(FLERR,arg[iarg+1]);
          if (cutoff <= 0.0)
              error->all(FLERR,"Illegal compute diamondlambda/atom command");
          cutsq = cutoff*cutoff;
          cutbig = cutoff*cutoff;
          iarg += 2;
      } else if (strcmp(arg[iarg],"cutoff_big") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal compute diamondlambda/atom command");
          double cutoff_big = force->numeric(FLERR,arg[iarg+1]);
          if (cutoff_big <= 0.0)
              error->all(FLERR,"Illegal compute diamondlambda/atom command");
          cutbig = cutoff_big*cutoff_big;
          iarg += 2;
      } else if (strcmp(arg[iarg],"hydrogen") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal compute diamondlambda/atom command");
          hydrogenId = force->numeric(FLERR,arg[iarg+1]);
          iarg += 2;
      } else if (strcmp(arg[iarg],"oxygen") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal compute diamondlambda/atom command");
          oxygenId = force->numeric(FLERR,arg[iarg+1]);
          iarg += 2;
      } else if (strcmp(arg[iarg],"deviation") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal compute diamondlambda/atom command");
          double deviation = force->numeric(FLERR,arg[iarg+1]);
          if (deviation < 0.0)
              error->all(FLERR,"Illegal compute diamondlambda/atom command");
          hydroDev=deviation;
          iarg += 2;
      } else if (strcmp(arg[iarg],"orderParameterOnly") == 0) {
          computeNucleiId=false;
          iarg += 1;
      } else if (strcmp(arg[iarg],"nucleiBiggest") == 0) {         
          biggest=true;
          iarg += 1;          
      } else if (strcmp(arg[iarg],"self") == 0) {
          if (arg[iarg+1][0] == 'g') compareDirection[iCondition] = true;
          else if (arg[iarg+1][0] == 'l') compareDirection[iCondition] = false;
          else error->all(FLERR,"Illegal compute diamondlambda/atom command");

          threshold[iCondition]=force->numeric(FLERR,arg[iarg+2]);
          iCondition++;
          iarg+=3;
      }
      else if (strcmp(arg[iarg],"hardNeighbour") == 0) {
          hardNeighbourDistance=force->numeric(FLERR,arg[iarg+1]);
          hardNeighbourCount=force->numeric(FLERR,arg[iarg+2]);
          iarg+=3;
      } else error->all(FLERR,"Illegal compute diamondlambda/atom command");
  }

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  comm_forward=2*(ndegree*2+1);
  NearestNeighNumber=NULL;          
  hydrogenBondNeigh=NULL;
  hydrogenBondNeigh_big=NULL;   
  qlmarray=NULL;
  qnvector = NULL;
  isSolid = NULL;
  nucleiID = NULL;
  maxneigh = 0;
  distsqO = NULL;
  distsqH = NULL;
  nearestO = NULL;
  nearestH = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeDiamondLambdaAtom::~ComputeDiamondLambdaAtom()
{
  memory->destroy(isSolid);
  memory->destroy(nucleiID);
  delete[] compareDirection;
  delete[] threshold;
  memory->destroy(NearestNeighNumber);           
  memory->destroy(hydrogenBondNeigh);
  memory->destroy(hydrogenBondNeigh_big);
  memory->destroy(qlmarray);
  memory->destroy(qnvector);
  memory->destroy(distsqO);
  memory->destroy(distsqH);
  memory->destroy(nearestO);
  memory->destroy(nearestH);
}

/* ---------------------------------------------------------------------- */

void ComputeDiamondLambdaAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute diamondlambda/atom requires a pair style be defined");
  if (cutsq == 0.0) cutsq = force->pair->cutforce * force->pair->cutforce;
  else if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR, "Compute diamondlambda/atom cutoff is longer than pairwise cutoff");
  if (rsoft<0) {
      error->all(FLERR, "Compute diamondlambda/atom rsoft is negative");
  }

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"diamondlambda/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute diamondlambda/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeDiamondLambdaAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */


int ComputeDiamondLambdaAtom::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc) {
    int i,j,m;
    m=0;
    for (i=0;i<n;i++) {
        j=list[i];
        if (packQlm) {
            for (int k=0;k<2*(ndegree*2+1);k++) {
                buf[m++]=qlmarray[j][k];
            }
        }
        if (packSolid) {
            buf[m++] = isSolid[j];
        }
        if (packNuclei) {
            buf[m++] = nucleiID[j];
        }
    }
    return m;
}
void ComputeDiamondLambdaAtom::unpack_forward_comm(int n, int first, double *buf) {
    int i,m,last;
    m=0;
    last=first+n;
    for (i=first;i<last;i++) {
        if (packQlm) {
            for (int k=0;k<2*(ndegree*2+1);k++) {
                qlmarray[i][k]=buf[m++];
            }
        }
        if (packSolid) {
            const bool newIsSolid = buf[m++];
            isSolid[i] = isSolid[i] || newIsSolid;
        }
        if (packNuclei) {
            const double oldNucleiId = nucleiID[i];
            const double newNucleiId = buf[m++];
            if (oldNucleiId == 0) {
              nucleiID[i] = newNucleiId;
            } else if (newNucleiId == 0) {
              nucleiID[i] = oldNucleiId;
            } else {
              nucleiID[i] = MIN(nucleiID[i], newNucleiId);
            }
        }
    }
}

void ComputeDiamondLambdaAtom::compute_peratom()
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow lambda parameter array if necessary

  if (atom->nlocal + atom->nghost > nmax) {
    memory->destroy(qlmarray);
    memory->destroy(qnvector);
    memory->destroy(isSolid);
    memory->destroy(nucleiID);
    nmax = atom->nmax;
    memory->create(qlmarray,nmax,(2*ndegree+1)*2,"diamondlambda/atom:qlmarray");
    memory->create(qnvector,nmax,"diamondlambda/atom:qnvector");
    memory->create(isSolid,nmax,"diamondlambda/atom:isSolid");
    memory->create(nucleiID,nmax,"diamondlambda/atom:nucleiID");
    if (computeNucleiId) {
        vector_atom = nucleiID;
    }
    else {
        vector_atom = qnvector;
    }
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  tagint *tag = atom->tag;

  memory->destroy(hydrogenBondNeigh);
  memory->destroy(hydrogenBondNeigh_big);
  //memory->create(hydrogenBondNeigh,inum,nnn,"diamondlambda/atom:hydrogenBondNeigh");   
  memory->create(hydrogenBondNeigh,inum,20,"diamondlambda/atom:hydrogenBondNeigh");  
  memory->create(hydrogenBondNeigh_big,inum,40,"diamondlambda/atom:hydrogenBondNeigh_big");  
  memory->destroy(NearestNeighNumber);              
  memory->create(NearestNeighNumber,inum,"diamondlambda/atom:NearestNeighNumber");     
  
  // compute lambda parameter for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double **x = atom->x;
  int *mask = atom->mask;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    hydrogenBondNeigh[ii][0]=-1;
    hydrogenBondNeigh_big[ii][0]=-1;
    if (oxygenId>=0&&atom->type[i]!=oxygenId) {
        continue;
    }
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      
      // insure distsq and nearest arrays are long enough

      if (jnum > maxneigh) {
        memory->destroy(distsqO);
        memory->destroy(distsqH);
        memory->destroy(nearestO);
        memory->destroy(nearestH);
        maxneigh = jnum;
        memory->create(distsqO,maxneigh,"diamondlambda/atom:distsqO");
        memory->create(distsqH,maxneigh,"diamondlambda/atom:distsqH");
        memory->create(nearestO,maxneigh,"diamondlambda/atom:nearestO");
        memory->create(nearestH,maxneigh,"diamondlambda/atom:nearestH");
      }

      // loop over list of all neighbors within force cutoff
      // distsq[] = distance sq to each
      // nearest[] = atom indices of neighbors
      // neighbor_big
      int ncountO = 0, ncountH = 0, nHydrongenBoundNeigh=0;
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        if (mask[j] & groupbit) {                       
                        delx = xtmp - x[j][0];
                        dely = ytmp - x[j][1];
                        delz = ztmp - x[j][2];

                        rsq = delx*delx + dely*dely + delz*delz;
                        if (rsq < cutbig) {
       	 			hydrogenBondNeigh_big[ii][nHydrongenBoundNeigh]=j;
        			nHydrongenBoundNeigh++;
        			if (nHydrongenBoundNeigh<40) {                  
            				hydrogenBondNeigh_big[ii][nHydrongenBoundNeigh]=-1;
                                }
                        }
                }
      }

      ncountO = 0, ncountH = 0;
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        if (mask[j] & groupbit) {                       
			delx = xtmp - x[j][0];
			dely = ytmp - x[j][1];
			delz = ztmp - x[j][2];

			rsq = delx*delx + dely*dely + delz*delz;
			if (rsq < cutsq) {
				if (oxygenId<0||atom->type[j]==oxygenId) {
					distsqO[ncountO] = rsq;
					nearestO[ncountO++] = j;
				}
				if (oxygenId>=0&&hydrogenId>=0&&atom->type[j]==hydrogenId) {
					distsqH[ncountH] = rsq;
					nearestH[ncountH++] = j;
				}
			}
		}
      }

      // use only nearest nnn neighbors
      NearestNeighNumber[ii] = ncountO;           
      
      //if (nnn<ncountO) {
      //   select2(nnn,ncountO,distsqO,nearestO);
      //    ncountO = nnn;
      //}
      //if (nnn<ncountH) {
      //    select2(nnn,ncountH,distsqH,nearestH);
      //    ncountH = nnn;
      //}
      
      nHydrongenBoundNeigh=0;
      for (jj = 0; jj < ncountO; jj++) {
        j = nearestO[jj];
        j &= NEIGHMASK;

        if (!hydrogenBond(i,j,ncountH)) {
            continue;
        }
        hydrogenBondNeigh[ii][nHydrongenBoundNeigh]=j;
        nHydrongenBoundNeigh++;
        if (nHydrongenBoundNeigh<20) {                  
            hydrogenBondNeigh[ii][nHydrongenBoundNeigh]=-1;
        }
      }
    }
  }
  for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      double *qlm= qlmarray[i];

      for (int m=0;m<2*ndegree+1;m++) {
          qlm[m*2]=0;
          qlm[m*2+1]=0;
          double sWeight=0;

          for (jj = 0; jj < 20; jj++) {     
              j = hydrogenBondNeigh[ii][jj];
              if (j==-1) {
                  break;
              }
              delx = xtmp - x[j][0];
              dely = ytmp - x[j][1];
              delz = ztmp - x[j][2];
              double r=sqrt(delx*delx+dely*dely+delz*delz);
              double rinv = 1.0/r;
              double weight=smearing(r);
              delx*=rinv;
              dely*=rinv;
              delz*=rinv;
              add_qlm_complex(m-ndegree,weight,delx,dely,delz,qlm+m*2,qlm+m*2+1);
              sWeight+=weight;
          }
          if (sWeight>0) {
              qlm[m*2]/=sWeight;
              qlm[m*2+1]/=sWeight;
          }
      }
  }
  packQlm=true;
  packSolid=false;
  packNuclei=false;
  comm->forward_comm_compute(this);
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      

      // loop over list of all neighbors within force cutoff
      // distsq[] = distance sq to each
      // nearest[] = atom indices of neighbors

      double sWeight=0;
      double usum=0,vsum=0;
      for (jj = 0; jj < 20; jj++) {           
          j = hydrogenBondNeigh[ii][jj];
          if (j == -1) {
              break;
          }

          double delx=atom->x[j][0]-atom->x[i][0];
          double dely=atom->x[j][1]-atom->x[i][1];
          double delz=atom->x[j][2]-atom->x[i][2];
          double weight=smearing(sqrt(delx*delx+dely*dely+delz*delz));
          add_qn_complex(i, j, weight, &usum, &vsum);
          sWeight+=weight;
      }
      if (sWeight>0) {
          usum/=sWeight;
          vsum/=sWeight;
      }
      qnvector[i]=usum;
  }
  if (!computeNucleiId) {
      return ;
  }
  for (int i = atom->nmax - 1; i >= 0; i -= 1) {
    isSolid[i] = 0;
    nucleiID[i] = 0;
  }
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if ((mask[i] & groupbit)) {
      if  (checkSolid(ii)&&NearestNeighNumber[ii]==nnn) { 
        isSolid[i] = 1;
        nucleiID[i] = tag[i];
      }
    }
  }
  //gn
  // if (biggest) {
  //	  for (ii = 0; ii<inum; ii++) {
  //		  i = ilist[ii];
  //		  if (isSolid[i] == 1){
  //			  for (jj = 0; jj < 20; jj++) {     
  //				  j = hydrogenBondNeigh[ii][jj];
  //				  if (j==-1) {
  //					  break;
  //				  isSolid[j] = 1;
  //				  nucleiID[j] = tag[j];
  //				  }
  //			  }
  //		  }
  //	  }
  // }
			  
	  
  packQlm=false;
  packSolid=true;
  packNuclei=false;
  comm->forward_comm_compute(this);
  packQlm=false;
  packSolid=false;
  packNuclei=true;
  comm->forward_comm_compute(this);

  int change,done,anychange;

  while (1) {
      change = 0;
      while (1) {
          done = 1;
          for (ii = 0; ii < inum; ii++) {
              i = ilist[ii];
              if (!(mask[i] & groupbit)) continue;
              if (!isSolid[i]) {
                  continue;
              }
              

              xtmp = x[i][0];
              ytmp = x[i][1];
              ztmp = x[i][2];
              for (jj = 0; jj < 40; jj++) {     
              //for (jj = ii; jj < inum; jj++){
                  j = hydrogenBondNeigh_big[ii][jj];
                  //j = ilist[jj];
                  if (j==-1) {
                      break;
                  }
                  if (nucleiID[i] == nucleiID[j]) {
                      continue;
                  }
                  if (!isSolid[j]) {
                      continue;
                  }

                  delx = xtmp - x[j][0];
                  dely = ytmp - x[j][1];
                  delz = ztmp - x[j][2];
                  rsq = delx*delx + dely*dely + delz*delz;
                  if (rsq < cutbig) {
                      int iMin = MIN(nucleiID[i],nucleiID[j]);
                      nucleiID[i] = nucleiID[j] = iMin;
                      done = 0;
                  }
              }
          }
          if (!done) change = 1;
          if (done) break;
      }
      MPI_Allreduce(&change,&anychange,1,MPI_INT,MPI_MAX,world);
      if (!anychange) break;
      packQlm=false;
      packSolid=false;
      packNuclei=true;
      comm->forward_comm_compute(this);
  }
  if (biggest) {
      for (ii = 0; ii < inum; ii++) {
          i = ilist[ii];
          if (!(mask[i] & groupbit)) continue;
          if (isSolid[i]) {
              continue;
          }
          nucleiID[i]=0;

          xtmp = x[i][0];
          ytmp = x[i][1];
          ztmp = x[i][2];

          for (jj = 0; jj < 20; jj++) {              
              j = hydrogenBondNeigh[ii][jj];
              if (j==-1) {
                  break;
              }
              if (nucleiID[i] == nucleiID[j]) {
                continue;
              }
              if (!isSolid[j]) {
                  continue;
              }

              delx = xtmp - x[j][0];
              dely = ytmp - x[j][1];
              delz = ztmp - x[j][2];
              rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < cutsq) {
                  if (nucleiID[i]==0||nucleiID[i]>nucleiID[j]) {
                      double iMin = nucleiID[j] +0.5;
                      nucleiID[i] = iMin;
                  }
              }
          }
      }
      MPI_Barrier(world);
  }
}

void ComputeDiamondLambdaAtom::add_qlm_complex(int m,double frr,double x,double y,double z,double *u,double *v) {
    double f=frr*polar_prefactor(ndegree,m,z);
    double phi=atan2(y,x);
    phi*=m;
    *u += f*cos(phi);
    *v += f*sin(phi);
}

// calculate lambda parameter using std::complex::pow function

inline void ComputeDiamondLambdaAtom::add_qn_complex(int i,int j,double frr, double *u, double *v) {
    double x=0,y=0;
    double normI=0,normJ=0;
    int m;
    for (m=0;m<2*ndegree+1;m++) {
        double xI=qlmarray[i][m*2],yI=qlmarray[i][m*2+1];
        double xJ=qlmarray[j][m*2],yJ=-qlmarray[j][m*2+1];
        double xD=xI*xJ-yI*yJ;
        double yD=xI*yJ+xJ*yI;
        double nID=xI*xI+yI*yI;
        double nJD=xJ*xJ+yJ*yJ;
        x+=xD;
        y+=yD;
        normI+=nID;
        normJ+=nJD;
    }
    double normIJ=sqrt(normI*normJ);
    if (normIJ>0) {
        *u+=x/normIJ*frr;
        *v+=y/normIJ*frr;
    }
}

/* ----------------------------------------------------------------------
   select2 routine from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   sort auxiliary array at same time
------------------------------------------------------------------------- */

#define SWAP(a,b)   tmp = a; a = b; b = tmp;
#define ISWAP(a,b) itmp = a; a = b; b = itmp;

/* ---------------------------------------------------------------------- */

void ComputeDiamondLambdaAtom::select2(int k, int n, double *arr, int *iarr)
{
  int i,ir,j,l,mid,ia,itmp;
  double a,tmp;

  arr--;
  iarr--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
        SWAP(arr[l],arr[ir])
        ISWAP(iarr[l],iarr[ir])
      }
      return;
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1])
      ISWAP(iarr[mid],iarr[l+1])
      if (arr[l] > arr[ir]) {
        SWAP(arr[l],arr[ir])
        ISWAP(iarr[l],iarr[ir])
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir])
        ISWAP(iarr[l+1],iarr[ir])
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1])
        ISWAP(iarr[l],iarr[l+1])
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      ia = iarr[l+1];
      for (;;) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        SWAP(arr[i],arr[j])
        ISWAP(iarr[i],iarr[j])
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      iarr[l+1] = iarr[j];
      iarr[j] = ia;
      if (j >= k) ir = j-1;
      if (j <= k) l = i;
    }
  }
}


/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDiamondLambdaAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  bytes += maxneigh * sizeof(double); 
  bytes += maxneigh * sizeof(int); 

  return bytes;
}

/* ----------------------------------------------------------------------
   polar prefactor for spherical harmonic Y_l^m, where 
   Y_l^m (theta, phi) = prefactor(l, m, cos(theta)) * exp(i*m*phi)
------------------------------------------------------------------------- */

double ComputeDiamondLambdaAtom::
polar_prefactor(int l, int m, double costheta) {
  const int mabs = abs(m);

  double prefactor = 1.0;
  for (int i=l-mabs+1; i < l+mabs+1; ++i)
    prefactor *= static_cast<double>(i);

  prefactor = sqrt(static_cast<double>(2*l+1)/(MY_4PI*prefactor))
    * associated_legendre(l,mabs,costheta);

  if ((m < 0) && (m % 2)) prefactor = -prefactor;
  //if (mabs%2) prefactor=-prefactor;

  return prefactor;
}

/* ----------------------------------------------------------------------
   associated legendre polynomial
------------------------------------------------------------------------- */

double ComputeDiamondLambdaAtom::
associated_legendre(int l, int m, double x) {
  if (l < m) return 0.0;

  double p(1.0), pm1(0.0), pm2(0.0);

  if (m != 0) {
    const double sqx = sqrt(1.0-x*x);
    for (int i=1; i < m+1; ++i)
      p *= static_cast<double>(2*i-1) * sqx;
  }

  for (int i=m+1; i < l+1; ++i) {
    pm2 = pm1;
    pm1 = p;
    p = (static_cast<double>(2*i-1)*x*pm1
         - static_cast<double>(i+m-1)*pm2) / static_cast<double>(i-m);
  }

  return p;
}
double ComputeDiamondLambdaAtom::smearing(double r) const {
    if (rsoft==0) {
        return 1;
    }
    else {
        return 1/(1+pow(r/rsoft,8));
    }
}

inline double dis(double x1,double y1,double z1,double x2,double y2,double z2) {
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}
bool ComputeDiamondLambdaAtom::hydrogenBond(int i,int j,int ncountH) const {
    if (hydroDev<0||hydrogenId<0) {
        return true;
    }
    bool ret=false;
    double **p=atom->x;
    double x1=p[i][0];
    double y1=p[i][1];
    double z1=p[i][2];
    double x2=p[j][0];
    double y2=p[j][1];
    double z2=p[j][2];
    double dc=dis(x1,y1,z1,x2,y2,z2)+hydroDev;
    int kk,k;
    //printf("%d\n",ncountH);
    for (kk=0;kk<ncountH;kk++) {
        k=nearestH[kk];
        k&=NEIGHMASK;
        double x3=atom->x[k][0];
        double y3=atom->x[k][1];
        double z3=atom->x[k][2];
        double d1=dis(x1,y1,z1,x3,y3,z3);
        double d2=dis(x2,y2,z2,x3,y3,z3);
        if (d1+d2<=dc) {
            ret=true;
        }
    }
    return ret;
}
bool ComputeDiamondLambdaAtom::checkSolid(int ii) const {
    int i=list->ilist[ii];
    if (oxygenId>=0&&atom->type[i]!=oxygenId) {
        return false;
    }
    if (hardNeighbourDistance>0) {
        double disSqr=hardNeighbourDistance*hardNeighbourDistance;
        int neighbourCount=0;
        double **p=atom->x;
        double x1=p[i][0];
        double y1=p[i][1];
        double z1=p[i][2];
        for (int jj = 0; jj < nnn; jj++) {            
            int j = hydrogenBondNeigh[ii][jj];
            if (j==-1) {
                break;
            }
            double x2=p[j][0];
            double y2=p[j][1];
            double z2=p[j][2];
            double xd=x2-x1;
            double yd=y2-y1;
            double zd=z2-z1;
            double currentDisSqr=xd*xd+yd*yd+zd*zd;
            if (currentDisSqr<=disSqr) {
                neighbourCount++;
            }
        }
        if (neighbourCount<hardNeighbourCount) {
            return false;
        }
    }
    bool ret=true;
    int m;
    for (m=0;m<nConditions;m++) {
        bool cd=compareDirection[m];
        double td=threshold[m];
        if (cd) {
            ret=ret&&qnvector[i]>td;
        }
        else {
            ret=ret&&qnvector[i]<td;
        }
        if (ret==false) {
            break;
        }
    }
    return ret;
}


