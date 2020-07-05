#include <linalg/tridiag.h>

using namespace std;

template <class T> int TTriDiagMatrix<T>::Solve(TVector<T>& tUnknowns,
						TVector<T>& tRHS)
{
  TVector<T> tGamma(_N+1);
  if(tUnknowns.Length()<_N+1)
    throw XError("Invalid unknowns!");
  
  T tBeta=_tB[0];
  if(tBeta==0)
    throw XError("Error solving tri-diagonal matrix--1");
  tUnknowns[0]=tRHS[0]/tBeta;
  
  T* a=&_tA[1];
  T* b=&_tB[1];
  T* c=&_tC[0];
  
  for(int j=1; j<=_N; j++)
    {
      tGamma[j]=*c/tBeta;
      tBeta=*b-*a*tGamma[j];
      if(tBeta==0)
	throw XError("Error solving tri-diagonal matrix--2");
      tUnknowns[j]=(tRHS[j]-*a*tUnknowns[j-1])/tBeta;
      a++;
      b++;
      c++;
    }
  
  for(int j=_N-1; j>=0; j--)
    {
      tUnknowns[j]-=tGamma[j+1]*tUnknowns[j+1];
    }
  
  return 1;
}

template <class T> int TCyclicTriDiagMatrix<T>::Solve(TVector<T>& tUnknowns, TVector<T>& tRHS)
{
  if(tUnknowns.Length()!=_N+1)
    throw XError("Invalid unknowns!");

  //TVector<T> tBB(_N+1);
  TVector<T> tZ(_N+1);
  TVector<T> tU(_N+1);

  T tGamma=_tB[0]*(-1);
  _tB[0]-=tGamma;
  //for(int i=1; i<_N; i++) tBB[i]=_tB[i];
  _tB[_N]-=_tC[_N]*_tA[0]/tGamma;

  //TTriDiagMatrix<T> t1(_tA, tBB, _tC, _N+1);

  TTriDiagMatrix<T>::Solve(tUnknowns, tRHS);

  tU[0]=tGamma;
  for(int i=1; i<_N; i++) tU[i]=0;
  tU[_N]=_tC[_N];

  TTriDiagMatrix<T>::Solve(tZ, tU);

  T tFactor=(tUnknowns[0]+_tA[0]*tUnknowns[_N]/tGamma)
    /(1+tZ[0]+_tA[0]*tZ[_N]/tGamma);
  for(int i=0; i<=_N; i++) tUnknowns[i]-=tFactor*tZ[i];

  return 1;
}
