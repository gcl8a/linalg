//01/21/03: revamped everything into one template class
//01/06/03: ...
//04/30/02: added return max(delta) to multigrid solver
//04/30/02: added epsilon to multigrid solver
//04/19/02: wrote comments
/*
multigrid.h and multigrid.cpp contain my multi-grid solvers
*/

#ifndef __LINALG_MULTIGRID_H
#define __LINALG_MULTIGRID_H

#include <template/array.h>
#include <linalg/linsolve.h>

using namespace std;

//enum RELAXATION_SCHEME {REDBLACK=1, MSI};
enum RELAXATION_SCHEME {REDBLACK=1, MSI, JACOBI, ADI, LINE_RELAX};


template <class T> class TMultiGrid
{
 protected:
  int _levels;
  TIArray<T> tiMatrices;
  
 public:
  TMultiGrid(int m, int n, int levels)
    : _levels(levels), tiMatrices(levels)
    {
      int k=0;
      while(k<levels)
	{
	  tiMatrices[k]=new T(m, n);
	  m=m/2;
	  n=n/2;
	  k++;
	}	    
    }
  
  ~TMultiGrid(void)
    {
      tiMatrices.Destroy();
    }

 protected:
  int Restrict(matrix&, matrix&, int, int);
  int RestrictOGrid(matrix&, matrix&, int, int);
  int Interpolate(matrix&, matrix&, int, int);
  int InterpolateOGrid(matrix&, matrix&, int, int);

 public:
  T& operator [] (int k) {return *tiMatrices[k];}

  PRECISION StepMG(int, int, int);
  PRECISION StepMGOGrid(int, int, int, RELAXATION_SCHEME);

 public:
  int PreDecompose(void);

  PRECISION SolveMG(int nu, int gamma) {return StepMG(nu, gamma, 0);}
  PRECISION SolveMGOGrid(int nu, int gamma) 
    {return StepMGOGrid(nu, gamma, 0, REDBLACK);}
  PRECISION SolveMGOGridMSI(int nu, int gamma) 
    {
      StepMGOGrid(nu, gamma, 0, MSI);
      return tiMatrices[0]->del.Max();//tiMatrices[0]->del);
    }
};

template <class T> int TMultiGrid<T>::Restrict(matrix& res, matrix& rhs, int M, int N)
  /*restricts the residual held in res to a half-sized matrix rhs
    M and N are the dimensions of the coarse grid, which runs [0..M][0..N]
   */
{
  for(int i=1; i<M; i++)
    {
      //ii and jj refer to fine grid
      int ii=i*2; 
      
      //homogeneous boundaries
      rhs[i][0]=0;
      rhs[i][N]=0;

      for(int j=1; j<N; j++)
	{
	  int jj=j*2;

	  //fancy-schmancy nine point stencil
	  rhs[i][j]=1.0*(0.0625*res[ii-1][jj-1]
			 +0.125*res[ii-1][jj]
			 +0.0625*res[ii-1][jj+1]
			 
			 +0.125*res[ii][jj-1]
			 +0.25*res[ii][jj]
			 +0.125*res[ii][jj+1]
			 
			 +0.0625*res[ii+1][jj-1]
			 +0.125*res[ii+1][jj]
			 +0.0625*res[ii+1][jj+1]
			 );
	}
    }
  
  //homogeneous boundaries
  for(int j=0; j<=N; j++)
    {
      rhs[0][j]=0;
      rhs[M][j]=0;
    }
 
  return 1;
}

template <class T> int TMultiGrid<T>::RestrictOGrid(matrix& res, matrix& rhs, int M, int N)
  /*same as ::Restrict, but assumes an overlapping boundary in j
    (i.e., rhs[i][0]=rhs[i][N]
   */
{
  for(int i=1; i<M; i++)
    {
      //ii and jj refer to fine grid
      int ii=i*2;
      int iim=ii-1;
      int iip=ii+1;

      for(int j=0; j<N; j++)
	{
	  int jj=j*2;
	  int jjm=jj ? jj-1 : 2*N-1;
	  int jjp=jj+1;  //jj<2*N ? jj+1 : 1; jj is always < 2*N

	  //fancy-schmancy nine point stencil
	  rhs[i][j]=0.0625*(res[iim][jjm]+res[iim][jjp]+res[iip][jjm]+res[iip][jjp])
	    +0.125*(res[iim][jj]+res[ii][jjm]+res[ii][jjp]+res[iip][jj])
	    +0.25*res[ii][jj];
	}

      //periodic boundary
      rhs[i][N]=rhs[i][0];
    }
  
  for(int j=0; j<=N; j++)
    {
      rhs[0][j]=0;
      rhs[M][j]=0;
    }
 
  return 1;
}

template <class T> int TMultiGrid<T>::Interpolate(matrix& cs, matrix& s, int M, int N)
  /*Interpolates the solution from a course grid, cs to a finer one, s
    M and N are the size of the coarse grid
   */
{
  for(int i=0; i<=M; i++)
    {
      //ii and jj refer to find grid
      int ii=i*2; 

      s[ii][0]+=1.0*cs[i][0];
      s[ii][2*N]+=1.0*cs[i][N];

      for(int j=0; j<=N; j++)
	{
	  int jj=j*2;

	  s[ii-1][jj-1]+=0.25*cs[i][j];
	  s[ii+1][jj-1]+=0.25*cs[i][j];
	  s[ii-1][jj+1]+=0.25*cs[i][j];
	  s[ii+1][jj+1]+=0.25*cs[i][j];
	  s[ii-1][jj]+=0.5*cs[i][j];
	  s[ii+1][jj]+=0.5*cs[i][j];
	  s[ii][jj-1]+=0.5*cs[i][j];
	  s[ii][jj+1]+=0.5*cs[i][j];
	  s[ii][jj]+=1.0*cs[i][j];
	}
    }

  return 1;
}

template <class T> int TMultiGrid<T>::InterpolateOGrid(matrix& cs, matrix& s, int M, int N)
  /*Same as ::Interpolate but assumes periodic in j
   */
{
  for(int i=1; i<M; i++)
    {
      //ii and jj refer to fine grid
      int ii=i*2;
      int iim=ii-1;
      int iip=ii+1;

      for(int j=0; j<N; j++)
	{
	  int jj=j*2;
	  int jjm=jj ? jj-1 : 2*N-1;
	  int jjp=jj+1; //jj is always < 2*N jj<2*N ? jj+1 : 1;

	  s[iim][jjm]+=0.25*cs[i][j];
	  s[iim][jj]+=0.5*cs[i][j];
	  s[iim][jjp]+=0.25*cs[i][j];
	  s[ii][jjm]+=0.5*cs[i][j];
	  s[ii][jj]+=1.0*cs[i][j];
	  s[ii][jjp]+=0.5*cs[i][j];
	  s[iip][jjm]+=0.25*cs[i][j];
	  s[iip][jj]+=0.5*cs[i][j];
	  s[iip][jjp]+=0.25*cs[i][j];
	}
    }

  //i=0 and i=M
  int MM=2*M;
  int MMm=MM-1;
  for(int j=0; j<N; j++)
    {
      int jj=j*2;
      int jjm=jj ? jj-1 : 2*N-1;
      int jjp=jj+1;

      s[0][jjm]+=0.5*cs[0][j];
      s[0][jjp]+=0.5*cs[0][j];
      s[0][jj]+=1.0*cs[0][j];
      s[1][jjm]+=0.25*cs[0][j];
      s[1][jjp]+=0.25*cs[0][j];
      s[1][jj]+=0.5*cs[0][j];

      s[MMm][jjm]+=0.25*cs[M][j];
      s[MMm][jjp]+=0.25*cs[M][j];
      s[MMm][jj]+=0.5*cs[M][j];
      s[MM][jjm]+=0.5*cs[M][j];
      s[MM][jjp]+=0.5*cs[M][j];
      s[MM][jj]+=1.0*cs[M][j];
    }

  for(int i=0; i<MM; i++)
    {
      s[i][2*N]=s[i][0];
      s[i][2*N]=s[i][0];
      s[i][2*N]=s[i][0];
    }
  
  return 1;
}

template <class T> int TMultiGrid<T>::PreDecompose(void)
{
  //pre-decomposes MSI solvers for use in multi-grid
  for(int k=0; k<_levels; k++)
    {
      tiMatrices[k]->PreDecompose();
    }
  return 1;
}

template <class T> PRECISION TMultiGrid<T>::StepMG(int nu, int gamma, int k)
{
  if(k<_levels-1)
    {
      for(int t=0; t<gamma; t++)
	{
	  //relax on this grid
	  tiMatrices[k]->RelaxRedBlack(nu, 0);

	  //restrict the residual
	  tiMatrices[k]->CalcResidual();

	  Restrict(tiMatrices[k]->res, tiMatrices[k+1]->rhs,
		   tiMatrices[k+1]->M, tiMatrices[k+1]->N);

	  tiMatrices[k+1]->sol.Zero();

	  //now call the solver for the next coarser grid
	  StepMG(nu, gamma, k+1);
      
	  Interpolate(tiMatrices[k+1]->sol, tiMatrices[k]->sol,
		   tiMatrices[k+1]->M, tiMatrices[k+1]->N);
	}
      
      //relax some more
      return tiMatrices[k]->RelaxRedBlack(nu, 0);
    }
  
  //if coarsest grid solve "directly"
  else return tiMatrices[k]->RelaxRedBlack(2*nu, 0);
}

template <class T> PRECISION TMultiGrid<T>::StepMGOGrid(int nu, int gamma, int k, RELAXATION_SCHEME scheme)
{
  if(k<_levels-1)
    {
      for(int t=0; t<gamma; t++)
	{
	  //relax on this grid
	  if(scheme==REDBLACK)
	    tiMatrices[k]->RelaxRedBlackOGrid(nu, 0);
	  else
	    tiMatrices[k]->SolveMSIOGrid(nu, 0);
	  
	  //restrict the residual
	  tiMatrices[k]->CalcResidualOGrid();

	  RestrictOGrid(tiMatrices[k]->res, tiMatrices[k+1]->rhs,
			tiMatrices[k+1]->M, tiMatrices[k+1]->N);
	  
	  tiMatrices[k+1]->sol.Zero();
	  
	  //now call the solver for the next coarser grid
	  StepMGOGrid(nu, gamma, k+1, scheme);
	  
	  InterpolateOGrid(tiMatrices[k+1]->sol, tiMatrices[k]->sol,
			   tiMatrices[k+1]->M, tiMatrices[k+1]->N);
	}
      
      //relax some more
      if(scheme==REDBLACK)
	return tiMatrices[k]->RelaxRedBlackOGrid(nu, 0);
      else
	return tiMatrices[k]->SolveMSIOGrid(nu, 0);
    }
  
  //if coarsest grid solve "directly"
  else if(scheme==REDBLACK)
    return tiMatrices[k]->RelaxRedBlackOGrid(nu, 0);
  else
    return tiMatrices[k]->SolveMSIOGrid(nu, 0);
}

template <class T> class TMultiGrid3D
{
 protected:
  int _levels;
  TIArray<T> tiMatrices;
  
 public:
  TMultiGrid3D(int I, int J, int K, int levels)
    : _levels(levels), tiMatrices(levels)
    {
      int l=0;
      while(l<levels)
	{
	  tiMatrices[l]=new T(I, J, K);
	  I=I/2;
	  J=J/2;
	  K=K/2;
	  l++;
	}	    
    }
  
  ~TMultiGrid3D(void)
    {
      tiMatrices.Destroy();
    }

 protected:
  int Restrict(matrix&, matrix&, int, int);
  int RestrictOGrid(matrix&, matrix&, int, int);
  int Interpolate(matrix&, matrix&, int, int);
  int InterpolateOGrid(matrix&, matrix&, int, int);

 public:
  T& operator [] (int l) {return *tiMatrices[l];}

  PRECISION StepMG(int, int, int);
  PRECISION StepMGOGrid(int, int, int, RELAXATION_SCHEME);

 public:
  int PreDecompose(void);

  //PRECISION SolveMG(int nu, int gamma) {return StepMG(nu, gamma, 0);}
  PRECISION SolveMGOGrid(int nu, int gamma, RELAXATION_SCHEME scheme)
    {return StepMGOGrid(nu, gamma, 0, scheme);}
};

#endif

template <class T> PRECISION TMultiGrid3D<T>::StepMGOGrid(int nu, int gamma, int k, RELAXATION_SCHEME scheme)
{
  if(k<_levels-1)
    {
      for(int t=0; t<gamma; t++)
	{
	  //relax on this grid
	  if(scheme==JACOBI)
	    tiMatrices[k]->RelaxGauss(nu);
	  else
	  {}

/*	  //restrict the residual
	  tiMatrices[k]->CalcResidualOGrid();

	  RestrictOGrid(tiMatrices[k]->res, tiMatrices[k+1]->rhs,
			tiMatrices[k+1]->M, tiMatrices[k+1]->N);
	  
	  tiMatrices[k+1]->sol.Zero();
	  
	  //now call the solver for the next coarser grid
	  StepMGOGrid(nu, gamma, k+1, scheme);
	  
	  InterpolateOGrid(tiMatrices[k+1]->sol, tiMatrices[k]->sol,
			   tiMatrices[k+1]->M, tiMatrices[k+1]->N);
*/	}
/*      
      //relax some more
      if(scheme==REDBLACK)
	return tiMatrices[k]->RelaxRedBlackOGrid(nu, 0);
      else
	return tiMatrices[k]->SolveMSIOGrid(nu, 0);
*/    }
  
  //if coarsest grid solve "directly"
  else if(scheme==JACOBI)
    return tiMatrices[k]->RelaxGauss(nu);
  else
    return 0;
}

