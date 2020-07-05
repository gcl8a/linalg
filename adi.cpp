#ifndef __LINALG_ADI_H
#include <linalg/adi.h>
#endif

TADIMatrix::TADIMatrix(int m, int n)
  :  M(m), N(n),
     a2(M+3,N+3,-1,-1), a4(M+3,N+3,-1,-1), a5(M+3,N+3,-1,-1), 
     a5ij(M+3,N+3,-1,-1), a5i(M+3,N+3,-1,-1), a5j(M+3,N+3,-1,-1),
     a6(M+3,N+3,-1,-1), a8(M+3,N+3,-1,-1),
     rhs0(M+3,N+3,-1,-1), rhs1(M+3,N+3,-1,-1),
     curr(M+3,N+3,-1,-1), next(M+3,N+3,-1,-1)
{}

int TADIMatrix::SolveADI(int kMax, double epsilon)
  /*
    Solves using an ADI routine.

    kMax is the maximum number of iterations. epsilon is the convergence limit
   */
{
  int t;

  for(t=0; t<kMax; t++)
    {
      //solve implicitly in j first
      matrix half(M+3, N+3, -1, -1);
      for(int i=0; i<=M; i++)
	{
	  TTriDiagMatrix<double> t(N+1);
	  TVector<double> tNewCol(N+1);
	  TVector<double> tRHS(N+1);
	  
	  for(int j=0; j<=N; j++)
            {
	      t._tA[j]=a2[i][j]/2;
	      t._tB[j]=a5ij[i][j]+a5[i][j]+a5j[i][j]/2;
	      t._tC[j]=a8[i][j]/2;
	      
	      tRHS[j]=(a5ij[i][j]-a5i[i][j]/2)*curr[i][j]
		-a4[i][j]*curr[i-1][j]/2
		-a6[i][j]*curr[i+1][j]/2
		+rhs0[i][j];
            }
	  
	  t.Solve(tNewCol, tRHS);
	  
	  for(int j=0; j<=N; j++)
            {
	      half[i][j]=tNewCol[j];
            }
	}
      
      //solve implicitly in i next, using n+1/2 values just calculated
      for(int j=0; j<=N; j++)
	{
	  TTriDiagMatrix<double> t(M+1);
	  TVector<double> tNewRow(M+1);
	  TVector<double> tRHS(M+1);
	  
	  for(int i=0; i<=M; i++)
            {
	      t._tA[i]=a4[i][j]/2;
	      t._tB[i]=a5ij[i][j]+a5[i][j]+a5i[i][j]/2;
	      t._tC[i]=a6[i][j]/2;
	      
	      tRHS[i]=(a5ij[i][j]-a5j[i][j]/2)*half[i][j]
		-a2[i][j]*half[i][j-1]/2
		-a8[i][j]*half[i][j+1]/2
		+rhs1[i][j];
            }
	  
	  t.Solve(tNewRow, tRHS);
	  
	  //copy rows to the new matrix
	  for(int i=0; i<=M; i++)
            {
	      next[i][j]=tNewRow[i];
            }
	}
      
      //find maximum delta
      double delta=0;
      for(int i=0; i<=M; i++)
      	{
	  for(int j=0; j<=N; j++)
	    {
	      double test=fabs(next[i][j]-curr[i][j]);
	      if(test>delta) delta=test;
	      curr[i][j]=next[i][j];
            }
	}
      
      if(delta<epsilon) break;
    }
  
  return t;
}

int TADIMatrix::SolveADIOGrid(int kMax, double epsilon)
  /*
    Solves using an ADI routine. The grid is assumed to be periodic in j

    kMax is the maximum number of iterations. epsilon is the convergence limit
   */
{
  int t;

  for(t=0; t<kMax; t++)
    {
      matrix half(M+3, N+3, -1, -1);

      //solve implicitly in j first
      for(int i=0; i<=M; i++)
	{
	  //n.b. 0 <=> N in the grid, but the solver assumes _no overlap_
	  TCyclicTriDiagMatrix<double> t(N);
	  TVector<double> tNewCol(N);
	  TVector<double> tRHS(N);
	  
	  for(int j=0; j<=N-1; j++)
            {
	      t._tA[j]=a2[i][j]/2;
	      t._tB[j]=a5ij[i][j]+a5[i][j]+a5j[i][j]/2;
	      t._tC[j]=a8[i][j]/2;
	      
	      tRHS[j]=(a5ij[i][j]-a5i[i][j]/2)*curr[i][j]
		-a4[i][j]*curr[i-1][j]/2
		-a6[i][j]*curr[i+1][j]/2
		+rhs0[i][j];
            }

	  t.Solve(tNewCol, tRHS);
	  
	  for(int j=0; j<N; j++)
            {
	      half[i][j]=tNewCol[j];
            }

	  //overlap
	  half[i][N]=half[i][0];
	}
      
      //solve implicitly in i next, using n+1/2 values just calculated
      for(int j=1; j<N; j++)
	{
	  TTriDiagMatrix<double> t(M+1);
	  TVector<double> tNewRow(M+1);
	  TVector<double> tRHS(M+1);
	  
	  for(int i=0; i<=M; i++)
            {
	      t._tA[i]=a4[i][j]/2;
	      t._tB[i]=a5ij[i][j]+a5[i][j]+a5i[i][j]/2;
	      t._tC[i]=a6[i][j]/2;
	      
	      tRHS[i]=(a5ij[i][j]-a5j[i][j]/2)*half[i][j]
		-a2[i][j]*half[i][j-1]/2
		-a8[i][j]*half[i][j+1]/2
		+rhs1[i][j];
            }
	  
	  t.Solve(tNewRow, tRHS);
	  
	  //copy rows to the new matrix
	  for(int i=0; i<=M; i++)
            {
	      next[i][j]=tNewRow[i];
            }
	}
      
      //j=0
      TTriDiagMatrix<double> t(M+1);
      TVector<double> tNewRow(M+1);
      TVector<double> tRHS(M+1);
      
      for(int i=0; i<=M; i++)
	{
	  t._tA[i]=a4[i][0]/2;
	  t._tB[i]=a5ij[i][0]+a5[i][0]+a5i[i][0]/2;
	  t._tC[i]=a6[i][0]/2;
	  
	  tRHS[i]=(a5ij[i][0]-a5j[i][0]/2)*half[i][0]
	    -a2[i][0]*half[i][N-1]/2
	    -a8[i][0]*half[i][1]/2
	    +rhs1[i][0];
	}
      
      t.Solve(tNewRow, tRHS);
      
      //copy rows to the new matrix
      for(int i=0; i<=M; i++)
	{
	  next[i][0]=tNewRow[i];
	  next[i][N]=next[i][0];
	}
      
      //find maximum delta
      double delta=0;
      for(int i=0; i<=M; i++)
      	{
	  for(int j=0; j<=N; j++)
	    {
	      double test=fabs(next[i][j]-curr[i][j]);
	      if(test>delta) delta=test;
	      curr[i][j]=next[i][j];
            }
	}
      
      if(delta<epsilon) break;
    }
  
  return t;
}

