#ifndef __VORTEX_H
#define __VORTEX_H

#ifndef __TEMPLATE_GRID_H
#include <template\grid.h>
#endif

#ifndef __TEMPLATE_MATRIX_H
#include <template\matrix.h>
#endif

#define PREVSTEP 0
#define THISSTEP 1
#define HALFSTEP 2
#define NEXTSTEP 3

/*vortex.h and vortex.cpp contain routines for solving the vorticity transport
equation for 2-D flow.

This equation has the form:

dW/dt = u(dW/dx) + v(dW/dy) + (1/Re) * Laplacian(W)

where:
	W=vorticity
	u=x velocity
	v=y velocity
	Re=Reynold's number
	Laplacian=d^2/dx^2 + d^2/dy^2 in cartesian coordinates
*/

enum MBoundaryType {DIRICHLET=1, NEUMANN, ROBBIN};

enum MFlowBoundaryType {NONE=0, NO_SLIP, MOVING_NO_SLIP};

class CFlowPoint
/*holds the data necessary to desribe a 2-D viscous flow using the vorticity
transport equation
*/
	{
	protected:
		CFlowPoint* _pcPointLeft;
		CFlowPoint* _pcPointRight;
		CFlowPoint* _pcPointTop;
		CFlowPoint* _pcPointBottom;

		//boundary type
		MFlowBoundaryType _mFlowBoundaryType;
		double _ffBoundaryValue;

		static double __ffGridDelta;
	public:
		MFlowBoundaryType GetFlowBoundaryType(void) {return _mFlowBoundaryType;}


		//velocities
		LTArray<double> _tU;
		LTArray<double> _tV;

		//streamfunction
		LTArray<double> _tStreamFunction;

		//vorticity
		LTArray<double> _tVorticity;

		CFlowPoint(void);

		double dVortdx(int);
		double dVortdy(int);
		double udVortdx(int nStep) {return _tU[nStep]*dVortdx;}
		double vdVortdy(int nStep) {return _tV[nStep]*dVortdy;}
		double d2Vortdx2(int);
		double d2Vortdy2(int);

		friend ostream& operator << (ostream&, const CFlowPoint&);
		friend class CViscousDomain;
	};

class CViscousDomain
/*class CDomain
*/
	{
	protected:
		const double _ffViscosity;
		const double _ffFreeStreamVelocity;
		const double _ffLength;
		const double _ffReynoldsNumber;

		const double _ffGridDelta;

		const double _ffBeta;
		const double _ffGamma;

		//storage for flow data -- four steps are needed for Adams-Bashforth
		LTGrid<CFlowPoint> _tFlowGrid;

		int SolveVorticityXDirection(int, int, int, double,
			MBoundaryType, MBoundaryType);
		int SolveVorticityYDirection(int, int, int, double,
			MBoundaryType, MBoundaryType);

		int SolveVorticityXDirectionADI(int, int, int, double, int);
		int SolveVorticityYDirectionADI(int, int, int, double, int);

		int CalculateVelocities(int);
	public:
		CViscousDomain(double, double, float, int, int);
		void Initialize(void);

		int NextStep(void);
		int SolveVorticityAB(double);
		int SolveVorticityADI(double);
		int SolveStreamFunctionGS(int, double);

		void SetFlowBoundary(int, int, int, int, MFlowBoundaryType, double=0);

		ostream& OutputU(ostream&, int) const;
		ostream& OutputV(ostream&, int) const;
		ostream& OutputStreamFunction(ostream&, int) const;
		ostream& OutputVorticity(ostream&, int) const;
		friend ostream& operator << (ostream&, const CViscousDomain&);
	};

#endif
