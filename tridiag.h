#ifndef __LINALG_TRIDIAG_H
#define __LINALG_TRIDIAG_H

#ifndef __TEMPLATE_XERROR_H
#include <template/xerror.h>
#endif

#ifndef __LINALG_VECTOR_H
#include <linalg/vector.h>
#endif

#ifndef __LINALG_MATRIX_H
#include <linalg/matrix.h>
#endif

using namespace std;

template <class T> class TTriDiagMatrix
/*implements a tri-diagonal matrix of the form:

| b0 c0 0  0  0  0  |
| a1 b1 c1 0  0  0  |
| 0  a2 b2 c2 0  0  |
...
| 0  0  0  0  aN bN |

N.B. that a1 and cN are not used
*/
{
 public:
  T* _tA;
  T* _tB;
  T* _tC;
  
 protected:
  const int _N;
  
 public:
  TTriDiagMatrix(void) : _N(-1)
    {}
  
  TTriDiagMatrix(int n)
    //tri-diagonal matrix of size n X n
    : _tA(new T[n]), _tB(new T[n]), _tC(new T[n]), _N(n-1)
    {}
  
  TTriDiagMatrix(const TTriDiagMatrix&); 				//not implemented
  TTriDiagMatrix& operator = (const TTriDiagMatrix&); 	//not implemented
  
/*    TTriDiagMatrix(TVector<T>& tA, TVector<T>& tB, TVector<T>& tC, int n) */
/*      : _tA(tA), _tB(tB), _tC(tC), _N(n-1) */
/*      {}  */
  
  virtual ~TTriDiagMatrix(void)
    {
      delete[] _tA;
      delete[] _tB;
      delete[] _tC;
    }

  virtual int Solve(TVector<T>&, TVector<T>&);
  
  //		friend class TMatrix<T>;
};

template <class T> class TCyclicTriDiagMatrix : public TTriDiagMatrix<T>
/*implements a cyclic tri-diagonal matrix of the form:
  
| b0 c0 0  0  0  a0 |
| a1 b1 c1 0  0  0  |
| 0  a2 b2 c2 0  0  |
...
| cN 0  0  0  aN bN |

N.B. that _N is the length of the array.
*/
{
 protected:
 public:
  TCyclicTriDiagMatrix(void) : TTriDiagMatrix<T>() {}
  TCyclicTriDiagMatrix(int n) : TTriDiagMatrix<T>(n) {}
  TCyclicTriDiagMatrix(TVector<T>& tA, TVector<T>& tB, TVector<T>& tC, int n)
    : TTriDiagMatrix<T>(tA, tB, tC, n) {}
  
  ~TCyclicTriDiagMatrix(void) {}
  virtual int Solve(TVector<T>&, TVector<T>&);
};


/*
template <class T> int TBandedTriDiagMatrix<T>::Solve(TVector<T>& tUnknowns,
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

  T tGamma=-_tB[0];
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
*/




template <class T> int TTriDiagMatrix<T>::Solve(TVector<T>& tUnknowns,
						TVector<T>& tRHS)
{
  TVector<T> tGamma(_N+1);
  if(tUnknowns.Length()<_N+1)
    throw XError("Invalid unknowns!");
  
  T tBeta=_tB[0];
  if(tBeta==0.0)
    throw XError("Error solving tri-diagonal matrix--1");
  
  tUnknowns[0]=(tRHS[0])/tBeta;
  
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

void makeEye(double& eye) {eye=1.0;}
void makeEye(dmatrix& eye)
{
	int size=eye.CountRows();
	eye=dmatrix::Eye(size); //dangerous, since we haven't test if eye is square
}

template <class T> int TCyclicTriDiagMatrix<T>::Solve(TVector<T>& tUnknowns, TVector<T>& tRHS)
{
	/////////////////////////////////////////////
	/////NOTE: B is altered by this function/////
	/////////////////////////////////////////////

  if(tUnknowns.Length()!=this->_N+1)
    throw XError("Invalid unknowns!");

  T eye(this->_tA[0]); //this sizes eye properly
  makeEye(eye);  //and this makes it equal to I (nb, if T is a double, we just get a scalar)

  //TVector<T> tBB(_N+1);
  TVector<T> tZ(this->_N+1);
  TVector<T> tU(this->_N+1);

  T tGamma=this->_tB[0]*(-1.0);
  this->_tB[0]-=tGamma; //the minus sign is ok...we're not trying to zero it
  //for(int i=1; i<_N; i++) tBB[i]=_tB[i];
  this->_tB[this->_N]-=this->_tC[this->_N]*(this->_tA[0]/tGamma);

  //TTriDiagMatrix<T> t1(_tA, tBB, _tC, _N+1);

  TTriDiagMatrix<T>::Solve(tUnknowns, tRHS);

  tU[0]=tGamma;
  for(int i=1; i<this->_N; i++) tU[i]=eye*0.0;
  tU[this->_N]=this->_tC[this->_N];

  TTriDiagMatrix<T>::Solve(tZ, tU);

  T tFactor=(tUnknowns[0]+(this->_tA[0]/tGamma)*tUnknowns[this->_N])
    /(eye+tZ[0]+(this->_tA[0]/tGamma)*tZ[this->_N]);

  for(int i=0; i<=this->_N; i++) tUnknowns[i]-=tZ[i]*tFactor;

  return 1;
}

class TBlockTriDiagMatrix
/*implements a tri-diagonal matrix of the form:

| b0 c0 0  0  0  0  |
| a1 b1 c1 0  0  0  |
| 0  a2 b2 c2 0  0  |
...
| 0  0  0  0  aN bN |

N.B. that a1 and cN are not used

where a,b,c are MxM matrices
*/
{
 public:
  dmatrix* _tA;
  dmatrix* _tB;
  dmatrix* _tC;
  
 protected:
	const int M; //size of blocks
	const int N; //number of blocks
  
 public:
//  TTriDiagMatrix(void) : N(-1)
//    {}
  
  TBlockTriDiagMatrix(int m, int n)
    : _tA(new dmatrix[n]), _tB(new dmatrix[n]), _tC(new dmatrix[n]), M(m), N(n-1)
    {}
  
  TBlockTriDiagMatrix(const TBlockTriDiagMatrix&); 					//not implemented
  TBlockTriDiagMatrix& operator = (const TBlockTriDiagMatrix&); 	//not implemented
  
/*    TTriDiagMatrix(TVector<T>& tA, TVector<T>& tB, TVector<T>& tC, int n) */
/*      : _tA(tA), _tB(tB), _tC(tC), _N(n-1) */
/*      {}  */
  
  virtual ~TBlockTriDiagMatrix(void)
    {
      delete[] _tA;
      delete[] _tB;
      delete[] _tC;
    }

//  virtual int Solve(TVector<T>&, TVector<T>&);
};



#endif

