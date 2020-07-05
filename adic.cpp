#ifndef __LINALG_ADIC_H
#include <linalg/adic.h>
#endif

TADIMatrix::TADIMatrix(int m, int n)
  :  _M(m+1), _N(n+1),
     a2(_M+1,_N+1), a4(_M+1,_N+1), a5(_M+1,_N+1), 
     a5ij(_M+1,_N+1), a5i(_M+1,_N+1), a5j(_M+1,_N+1),
     a6(_M+1,_N+1), a8(_M+1,_N+1),
     rhs0(_M+1,_N+1), rhs1(_M+1,_N+1),
     curr(_M+1,_N+1), next(_M+1,_N+1)
{}

int TADIMatrix::SolveCGridADI(int nTrailingEdge, double tSing0,
			      double tSing1, int nMaxIter, double ffConv)
/*
*/
{
  int t;
  for(int t=0; t<nMaxIter; t++)
    {
      TMatrix<double> tHalfStep(_M+1,_N+1);
      //solve implicitly in j first
      
      //special case at i=1 because of periodic nature
      int i=1;
      TTriDiagMatrix<double> t(nTrailingEdge);
      TVector<double> tNewCol(nTrailingEdge);
      TVector<double> tRHS(nTrailingEdge);
      
      for(int j=1; j<nTrailingEdge; j++)
	{
	  t._tA[j-1]=a2[i][j];
	  t._tB[j-1]=a5ij[i][j]+a5j[i][j]+a5[i][j];
	  t._tC[j-1]=a8[i][j];
	  
	  tRHS[j-1]=(a5ij[i][j]-a5i[i][j])*curr[i][j]
	    -a4[i][j]*curr[i+1][_N-j]
	    -a6[i][j]*curr[i+1][j]
	    +rhs0[i][j];
	}
      
      //last equation is for the singularity
      t._tA[nTrailingEdge-1]=0.0;
      t._tB[nTrailingEdge-1]=1.0;
      t._tC[nTrailingEdge-1]=0.0;
      tRHS[nTrailingEdge-1]=tSing0;
      
      t.Solve(tNewCol, tRHS);
      
      //copy result to _both_ boundaries at the wake, but not the trailing edge
      for(int j=1; j<nTrailingEdge; j++)
	{
	  tHalfStep[i][j]=tNewCol[j-1];
	  tHalfStep[i][_N-j]=tNewCol[j-1];
	}
      
      //don't bother solving for airfoil now-just copy the values because we have
      //Dirichlet bc's
      for(int j=nTrailingEdge; j<=_N-nTrailingEdge; j++)
	{
	  tHalfStep[i][j]=curr[i][j];
	}
      
      //solve for all of the other columns with no special favors
      for(i=2; i<_M; i++)
	{
	  TTriDiagMatrix<double> t(_N-1);
	  TVector<double> tNewCol(_N-1);
	  TVector<double> tRHS(_N-1);
	  
	  for(int j=1; j<_N; j++)
            {
	      t._tA[j-1]=a2[i][j];
	      t._tB[j-1]=a5ij[i][j]+a5j[i][j]+a5[i][j];
	      t._tC[j-1]=a8[i][j];
	      
	      tRHS[j-1]=(a5ij[i][j]-a5i[i][j])*curr[i][j]
		-a4[i][j]*curr[i-1][j]
		-a6[i][j]*curr[i+1][j]
		+rhs0[i][j];
            }
	  
	  t.Solve(tNewCol, tRHS);
	  
	  for(int j=1; j<_N; j++)
            {
	      tHalfStep[i][j]=tNewCol[j-1];
            }
	}
      
      //solve implicitly in i next, using n+1/2 values just calculated
      
      //in the downstream region, we solve one tri-diag for the whole shootin' match
      for(int j=1; j<nTrailingEdge; j++)
	{
	  TTriDiagMatrix<double> t(2*_M-3);
	  TVector<double> tNewRow(2*_M-3);
	  TVector<double> tRHS(2*_M-3);
	  
	  int i=1;
	  t._tA[_M-i-1]=a6[i][j];
	  t._tB[_M-i-1]=a5ij[i][j]+a5i[i][j]+a5[i][j];
	  t._tC[_M-i-1]=a4[i][j];  //could be a6[i][_N-j]
	  
	  //if at trailing edge-1, we must use the singularity value at the trailing edge
	  tRHS[_M-i-1]= (j==nTrailingEdge-1) ?
	    (a5ij[i][j]-a5j[i][j])*tHalfStep[i][j]
	    -a2[i][j]*tHalfStep[i][j-1]
	    -a8[i][j]*tSing1
	    +rhs1[i][j]
	    
	    : (a5ij[i][j]-a5j[i][j])*tHalfStep[i][j]
	    -a2[i][j]*tHalfStep[i][j-1]
	    -a8[i][j]*tHalfStep[i][j+1]
	    +rhs1[i][j];
	  
	  for(i=2; i<_M; i++)
            {
	      t._tA[_M-i-1]=a6[i][j];
	      t._tB[_M-i-1]=a5ij[i][j]+a5i[i][j]+a5[i][j];
	      t._tC[_M-i-1]=a4[i][j];
	      
	      tRHS[_M-i-1]=(a5ij[i][j]-a5j[i][j])*tHalfStep[i][j]
		-a2[i][j]*tHalfStep[i][j-1]
		-a8[i][j]*tHalfStep[i][j+1]
		+rhs1[i][j];
	      
	      t._tA[_M+i-3]=a4[i][_N-j];
	      t._tB[_M+i-3]=a5ij[i][_N-j]+a5i[i][_N-j]+a5[i][_N-j];
	      t._tC[_M+i-3]=a6[i][_N-j];
	      
	      tRHS[_M+i-3]=(a5ij[i][_N-j]-a5j[i][_N-j])*tHalfStep[i][_N-j]
		-a2[i][_N-j]*tHalfStep[i][_N-j-1]
		-a8[i][_N-j]*tHalfStep[i][_N-j+1]
		+rhs1[i][_N-j];
            }
	  
	  t.Solve(tNewRow, tRHS);
	  
	  //copy rows to the new matrix
	  for(i=1; i<_M; i++)
            {
	      next[i][j]=tNewRow[_M-i-1];
	      next[i][_N-j]=tNewRow[_M+i-3];
            }
	}
      
      //now do the part for the airfoil -- solve tri-diags
      for(int j=nTrailingEdge; j<=_N-nTrailingEdge; j++)
	{
	  TTriDiagMatrix<double> t(_M-1);
	  TVector<double> tNewRow(_M-1);
	  TVector<double> tRHS(_M-1);
	  
	  for(int i=1; i<_M; i++)
            {
	      t._tA[i-1]=a4[i][j];
	      t._tB[i-1]=a5ij[i][j]+a5i[i][j]+a5[i][j];
	      t._tC[i-1]=a6[i][j];
	      
	      tRHS[i-1]=(a5ij[i][j]-a5j[i][j])*tHalfStep[i][j]
		-a2[i][j]*tHalfStep[i][j-1]
		-a8[i][j]*tHalfStep[i][j+1]
		+rhs1[i][j];
            }
	  
	  t.Solve(tNewRow, tRHS);
	  
	  //copy rows to the new matrix
	  for(int i=1; i<_M; i++)
            {
	      next[i][j]=tNewRow[i-1];
            }
	}
      
      //copy cyclics and verify that i==0 is equal
      for(int j=1; j<nTrailingEdge; j++)
	{
	  if(next[1][_N-j]!=next[1][j])
            throw XError("ADI error in wake!");
	  
	  //		nextdata[0][j].vort=nextdata[2][_N-j].vort;
	  //		nextdata[0][_N-j].vort=nextdata[2][j].vort;
	}
      
      //ugly
      double conv=0;
      for(int i=1; i<_M; i++)
      	{
	  for(int j=1; j<_N; j++)
	    {
	      double test=fabs(curr[i][j]-next[i][j]);
	      if(test>conv) conv=test;
	      curr[i][j]=next[i][j];
            }
	}
      
      if(conv<ffConv) break;
    }
  
  return t;
}

int TADIMatrix::SolveCGridCN(int nTrailingEdge, double tSing0,
			     double tSing1, int nMaxIter, double ffConv)
  /*solves a c grid using Crank-Nicolson
*/
{
  TMatrix<double> tNext(_M+1,_N+1);

  int t;
  for(int t=0; t<nMaxIter; t++)
    {
      //special case at i=1 because of periodic nature
      int i=1;
      
      for(int j=2; j<nTrailingEdge-1; j++)
	{
	  tNext[i][j]=curr[i][j]
	    -a2[i][j]*(curr[i][j-1]+next[i][j-1])
	    -a4[i][j]*(curr[i+1][_N-j]+next[i+1][_N-j])
	    -a6[i][j]*(curr[i+1][j]+next[i+1][j])
	    -a8[i][j]*(curr[i][j+1]+next[i][j+1])
	    
	    -(a5i[i][j]+a5j[i][j])*(curr[i][j]+next[i][j]);
	  tNext[i][_N-j]=tNext[i][j];
	}
      
      int j=nTrailingEdge-1;
      tNext[i][j]=curr[i][j]
	-a2[i][j]*(curr[i][j-1]+next[i][j-1])
	-a4[i][j]*(curr[i+1][_N-j]+next[i+1][_N-j])
	-a6[i][j]*(curr[i+1][j]+next[i+1][j])
	-a8[i][j]*(tSing0+tSing1)
	    
	-(a5i[i][j]+a5j[i][j])*(curr[i][j]+next[i][j]);
      tNext[i][_N-j]=tNext[i][j];

      for(int j=nTrailingEdge; j<=_N-nTrailingEdge; j++)
	{
	  tNext[i][j]=rhs1[i][j];
	}

      for(int i=2; i<_M; i++)
	{
	  for(int j=2; j<_N-1; j++)
	    {
	      tNext[i][j]=curr[i][j]
		-a2[i][j]*(curr[i][j-1]+next[i][j-1])
		-a4[i][j]*(curr[i-1][j]+next[i-1][j])
		-a6[i][j]*(curr[i+1][j]+next[i+1][j])
		-a8[i][j]*(curr[i][j+1]+next[i][j+1])

		-(a5i[i][j]+a5j[i][j])*(curr[i][j]+next[i][j]);
	    }
	}
      
      for(int i=1; i<_M; i++)
	{
	  tNext[i][1]=rhs1[i][1];
	  tNext[i][_N-1]=rhs1[i][_N-1];
	}

      //copy cyclics and verify that i==0 is equal
      for(int j=1; j<nTrailingEdge; j++)
	{
	  if(tNext[1][_N-j]!=tNext[1][j])
            throw XError("C-N error in wake!");
	  
	  //		nextdata[0][j].vort=nextdata[2][_N-j].vort;
	  //		nextdata[0][_N-j].vort=nextdata[2][j].vort;
	}

      double conv=0;
      for(int i=1; i<_M; i++)
      	{
	  for(int j=1; j<_N; j++)
	    {
	      double test=fabs(tNext[i][j]-next[i][j]);
	      if(test>conv) conv=test;
	      next[i][j]=tNext[i][j];
            }
	}
      
      if(conv<ffConv) break;
    }
      
  return t;
}
