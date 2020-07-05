#ifndef __LINALG_ADIC_H
#define __LINALG_ADIC_H

#ifndef __LINALG_TRIDIAG_H
#include <linalg/tridiag.h>
#endif

#ifndef __LINALG_MATRIX_H
#include <linalg/matrix.h>
#endif

class TADIMatrix
{
 protected:
  const int _M, _N;
  
 public:
  TMatrix<double> a2;
  TMatrix<double> a4;
  TMatrix<double> a5;
  TMatrix<double> a5ij;
  TMatrix<double> a5i;
  TMatrix<double> a5j;
  TMatrix<double> a6;
  TMatrix<double> a8;
  
  TMatrix<double> rhs0;
  TMatrix<double> rhs1;
  TMatrix<double> curr;
  TMatrix<double> next;

 public:
  TADIMatrix(int m, int n);
  
  int SolveCyclicADI(double, int);
  int SolveCGridADI(int, double, double, int, double);
  int SolveCGridCN(int, double, double, int, double);
  int Solve(int);
};

#endif




