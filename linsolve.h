//01/21/03: maintenance -- made each diag matrix more similar
//11/15/02: removed #define ALPHA. alpha now varies for each size matrix
//11/15/02: edited for consistancy. changed deficit to residual
//11/14/2002: fixed minor return value warning
//05/06/02: corrected some errors in ESI routine, set ALPHA=0.5 by default
//04/30/02: changed /E to *E in MSI routines
//04/30/02: MSI solver (for T13) skips finding max delta if epsilon=0
//04/19/02: commented out some of the confusing multigrid stuff

/*
linsolve.h and linsolve.cpp provide numerous solvers for elliptic problems, such as the Poisson equation for stream function or pressure.
*/
#define PRECISION double

#ifndef __LINALG_LINSOLVE_H
#define __LINALG_LINSOLVE_H

#include <linalg/matrix.h>

using namespace std;

class TDiagMatrix
{
public:
	const int M, N;

	matrix sol; 
	matrix rhs; 
	matrix res; //res=rhs-A*sol
	matrix del; //delta solution

public:
	TDiagMatrix(int m, int n);
	~TDiagMatrix(void) {}

	//unimplemented -- copying could be bad
	TDiagMatrix(const TDiagMatrix&);
	TDiagMatrix& operator = (const TDiagMatrix&);
};

class TPentaDiagMatrix : public TDiagMatrix
{
public:
	matrix a2;
	matrix a4;
	matrix a5;
	matrix a6;
	matrix a8;

	matrix B;
	matrix C;
	matrix D;
	matrix E;
	matrix F;
	matrix G;
	matrix H;

public:
	TPentaDiagMatrix(int m, int n);
	~TPentaDiagMatrix(void) {}

	//unimplemented -- copying could be bad
	TPentaDiagMatrix(const TPentaDiagMatrix&);
	TPentaDiagMatrix& operator = (const TPentaDiagMatrix&);

	int SetupMSI(void);
	int PreDecompose(void) {return PartialLUDecompOGrid(0.9);}
	int PartialLUDecomp(float);
	int PartialLUDecompOGrid(float);
	int StepMSI(void);
	int SolveMSI(int, PRECISION);
	int SolveMSIOGrid(int, PRECISION);

	int StepGS(void);
	PRECISION RelaxRedBlack(int, PRECISION);
	PRECISION RelaxRedBlackOGrid(int, PRECISION);

	int CalcResidual(void);
	int CalcResidualOGrid(void);
	int Output(char*);
};

class TNonaDiagMatrix : public TPentaDiagMatrix
{
public:
	matrix a1;
	matrix a3;
	matrix a7;
	matrix a9;

	matrix A;
	matrix U;

public:
	TNonaDiagMatrix(int m, int n);
	~TNonaDiagMatrix(void) {}

	//unimplemented -- copying could be bad
	TNonaDiagMatrix(const TNonaDiagMatrix&);
	TNonaDiagMatrix& operator = (const TNonaDiagMatrix&);

	int SetupMSI(void);
	int PreDecompose(void) {return PartialLUDecompOGrid(0.9);}
	int PartialLUDecomp(float);
	int PartialLUDecompOGrid(float);
	int StepMSI(void);
	int SolveMSI(int, PRECISION);
	int SolveMSIOGrid(int, PRECISION);

	int StepGS(void);
	PRECISION RelaxRedBlack(int, PRECISION);
	PRECISION RelaxRedBlackOGrid(int, PRECISION);

	int CalcResidual(void);
	int CalcResidualOGrid(void);
	int Output(char*);
};

class T13DiagMatrix : public TNonaDiagMatrix
{
public:
	matrix ass;
	matrix aww;
	matrix aee;
	matrix ann;

	matrix P;
	matrix Q;
	matrix R;
	matrix S;
public:
	T13DiagMatrix(int m, int n);
	~T13DiagMatrix(void) {}

	//unimplemented -- copying could be bad
	T13DiagMatrix(const T13DiagMatrix&);
	T13DiagMatrix& operator = (const T13DiagMatrix&);

	int SetupMSI(void);
	int PreDecompose(void) {return PartialLUDecompOGrid(0.5);}
	int StepMSI(void);
	int PartialLUDecomp(float);
	int PartialLUDecompOGrid(float);
	int SolveMSI(int, PRECISION);
	PRECISION SolveMSIOGrid(int, PRECISION);

	int StepGS(void);
	PRECISION RelaxRedBlack(int, PRECISION);
	PRECISION RelaxRedBlackOGrid(int, PRECISION);

	int CalcResidual(void);
	int CalcResidualOGrid(void);
	int Output(char*);
};


class TBanded3D
	/*
	This is called TBanded3D for lack of a better term. It is for solving a 3D system
	of equations where the cells are related to each other by matrices.

	Hard coded for four variables and 4x4 matrix of 'interconnections', in the
	orthogonal direction only.

	The system of equations is:

	a0*q0 + aW*qW + aE*qE + ... = rhs

	where the rhs and q's are 4-vectors and the a's are 4x4 matrices
	*/
{
protected:
	int iMax, jMax, kMax;
public:
	T4DMatrix x;

	T5DMatrix a0, aW, aE, aS, aN, aU, aD;
	T4DMatrix rhs;

	//T4DMatrix R;  //equals R + Q_mn (a portion of S)
	//T4DMatrix f;

	TBanded3D(int I, int J, int K);
	int RelaxGauss(int);
	int LineRelax(void);
	int LineRelaxI(int, int);
	int LineRelaxJ(int, int);
	int LineRelaxK(int, int);

	int Zero(void);

/*	int aCalcRHS(void) 
	{
		for(int i=0; i<=iMax; i++)
			for(int j=0; j<=jMax; j++)
				for(int k=0; k<=kMax; k++)
					rhs[i][j][k];//;=f[i][j][k]+R[i][j][k]-R0[i][j][k];

		return 0;
	}*/
};

#endif
