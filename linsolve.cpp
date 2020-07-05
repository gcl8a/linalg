#ifndef __LINALG_LINSOLVE_H
#include <linalg/linsolve.h>
#endif

#include <fstream>

#include <linalg/vector.h>
#include <linalg/matrix.h>
#include <linalg/tridiag.h>

using namespace std;

TDiagMatrix::TDiagMatrix(int m, int n)
  : M(m), N(n),
    sol(m+5, n+5, -2, -2), rhs(m+5, n+5, -2, -2), res(m+5, n+5, -2, -2), 
    del(m+5, n+5, -2, -2)
{}

TPentaDiagMatrix::TPentaDiagMatrix(int m, int n)
  : TDiagMatrix(m, n),
    a2(m+5, n+5, -2, -2), a4(m+5, n+5, -2, -2), a5(m+5, n+5, -2, -2), 
    a6(m+5, n+5, -2, -2), a8(m+5, n+5, -2, -2)
{}

TNonaDiagMatrix::TNonaDiagMatrix(int m, int n)
  : TPentaDiagMatrix(m, n),
    a1(m+5, n+5, -2, -2), a3(m+5, n+5, -2, -2), 
    a7(m+5, n+5, -2, -2), a9(m+5, n+5, -2, -2)
{}

T13DiagMatrix::T13DiagMatrix(int m, int n)
  : TNonaDiagMatrix(m, n),
    ass(m+5, n+5, -2, -2), aww(m+5, n+5, -2, -2), 
    aee(m+5, n+5, -2, -2), ann(m+5, n+5, -2, -2)
{}

int TPentaDiagMatrix::PartialLUDecomp(float alpha)
{
  int il=0;
  int iu=M;
  int jl=0;
  int ju=N;

  for(int i=il; i<=iu; i++)
    {
      for(int j=jl; j<=ju; j++)
	{
	  B[i][j]=(a4[i][j])
	    /(1.0-alpha*F[i-1][j]*F[i-1][j+1]);
	  C[i][j]=-B[i][j]*F[i-1][j];
	  D[i][j]=(a2[i][j]-B[i][j]*G[i-1][j])
	    /(1.0+2.0*alpha*G[i][j-1]);
	  
	  PRECISION phi1=C[i][j]*F[i-1][j+1];
	  PRECISION phi4=D[i][j]*G[i][j-1];
	  
	  E[i][j]=a5[i][j]-B[i][j]*H[i-1][j]-C[i][j]*G[i-1][j+1]-D[i][j]*F[i][j-1]
	    +alpha*2.0*(phi1+phi4);
	  E[i][j]=1.0/E[i][j];
	  F[i][j]=(a8[i][j]-C[i][j]*H[i-1][j+1]-2.0*alpha*(phi1))*E[i][j];
	  G[i][j]=(-D[i][j]*H[i][j-1])*E[i][j];
	  H[i][j]=(a6[i][j]-alpha*phi4)*E[i][j];
	}
    }
  
  return 1;
}

int TNonaDiagMatrix::PartialLUDecomp(float alpha)
{
  int il=0;
  int iu=M;
  int jl=0;
  int ju=N;

  for(int i=il; i<=iu; i++)
    {
      for(int j=jl; j<=ju; j++)
	{
	  A[i][j]=a1[i][j];
	  B[i][j]=(a4[i][j]-A[i][j]*F[i-1][j-1]-alpha*a7[i][j]*F[i-1][j+1])
	    /(1.0-alpha*F[i-1][j]*F[i-1][j+1]);
	  C[i][j]=a7[i][j]-B[i][j]*F[i-1][j];
	  D[i][j]=(a2[i][j]-A[i][j]*H[i-1][j-1]-B[i][j]*G[i-1][j]
		   -2*alpha*A[i][j]*G[i-1][j-1])
	    /(1.0+2.0*alpha*G[i][j-1]);
	  
	  PRECISION phi1=C[i][j]*F[i-1][j+1];
	  PRECISION phi2=A[i][j]*G[i-1][j-1];
	  PRECISION phi3=C[i][j]*U[i-1][j+1];
	  PRECISION phi4=D[i][j]*G[i][j-1];
	  
	  E[i][j]=a5[i][j]-A[i][j]*U[i-1][j-1]-B[i][j]*H[i-1][j]
	    -C[i][j]*G[i-1][j+1]-D[i][j]*F[i][j-1]
	    +alpha*(2.0*phi1+phi2+phi3+2.0*phi4);
	  E[i][j]=1.0/E[i][j];
	  F[i][j]=(a8[i][j]-B[i][j]*U[i-1][j]-C[i][j]*H[i-1][j+1]
		   -2.0*alpha*(phi1+phi3))*E[i][j];
	  G[i][j]=(a3[i][j]-D[i][j]*H[i][j-1])*E[i][j];
	  H[i][j]=(a6[i][j]-D[i][j]*U[i][j-1]-alpha*phi4)*E[i][j];
	  U[i][j]=a9[i][j]*E[i][j];
	}
    }
  
  return 1;
}

int T13DiagMatrix::PartialLUDecomp(float alpha)
{
  int il=0;
  int iu=M;
  int jl=0;
  int ju=N;
  
  for(int i=il; i<=iu; i++)
    {
      for(int j=jl; j<=ju; j++)
	{
	  P[i][j]=aww[i][j];
	  
	  PRECISION phi1=P[i][j]*F[i-2][j];
	  PRECISION phi2=P[i][j]*R[i-2][j];
	  
	  A[i][j]=a1[i][j]-P[i][j]*G[i-2][j];
	  B[i][j]=(a4[i][j]-P[i][j]*H[i-2][j]-A[i][j]*F[i-1][j-1]-alpha*
		   (2*phi1+2*phi2+(F[i-1][j+1]+R[i-1][j+1])
		    *(a7[i][j]-P[i][j]*U[i-2][j]-A[i][j]*R[i-1][j-1])))
	    /(1.0+alpha*(R[i-1][j]-F[i-1][j]*(F[i-1][j+1]+R[i-1][j+1])));
	  C[i][j]=a7[i][j]-P[i][j]*U[i-2][j]-A[i][j]*R[i-1][j-1]-B[i][j]*F[i-1][j];
	  
	  PRECISION phi3=B[i][j]*R[i-1][j]+C[i][j]*F[i-1][j+1];
	  PRECISION phi4=C[i][j]*R[i-1][j+1];
	  
	  Q[i][j]=ass[i][j]-A[i][j]*G[i-1][j-1];
	  
	  PRECISION phi5=Q[i][j]*G[i][j-2];
	  PRECISION phi7=Q[i][j]*S[i][j-2];
	  
	  D[i][j]=(a2[i][j]-A[i][j]*H[i-1][j-1]-B[i][j]*G[i-1][j]-Q[i][j]*F[i][j-2]
		   -alpha*(3*phi5+2*phi7)-alpha*2*Q[i][j]*H[i][j-2])
	    /(1.0+alpha*(2*G[i][j-1]+S[i][j-1]));
	  
	  PRECISION phi6=Q[i][j]*H[i][j-2]+D[i][j]*G[i][j-1];
	  PRECISION phi8=D[i][j]*S[i][j-1];
	  
	  E[i][j]=a5[i][j]-P[i][j]*S[i-2][j]-A[i][j]*U[i-1][j-1]-B[i][j]
	    *H[i-1][j]-C[i][j]*G[i-1][j+1]
	    -Q[i][j]*R[i][j-2]-D[i][j]*F[i][j-1]+alpha
	    *(2*phi1+3*phi2+2*phi3+3*phi4+3*phi5+2*phi6+3*phi7+2*phi8);
	  E[i][j]=1.0/E[i][j];
	  F[i][j]=(a8[i][j]-B[i][j]*U[i-1][j]-C[i][j]*H[i-1][j+1]-D[i][j]
		   *R[i][j-1]-alpha*(phi1+2*phi2+2*phi3+3*phi4))*E[i][j];
	  R[i][j]=(ann[i][j]-C[i][j]*U[i-1][j+1])*E[i][j];
	  G[i][j]=(a3[i][j]-A[i][j]*S[i-1][j-1]-Q[i][j]*U[i][j-2]-D[i][j]
		   *H[i][j-1])*E[i][j];
	  H[i][j]=(a6[i][j]-B[i][j]*S[i-1][j]-D[i][j]*U[i][j-1]
		   -alpha*(phi5+phi6+2*phi7+2*phi8))*E[i][j];
	  U[i][j]=(a9[i][j]-C[i][j]*S[i-1][j+1])*E[i][j];
	  S[i][j]=aee[i][j]*E[i][j];
	}
    }
  
  return 1;
}

int TPentaDiagMatrix::SolveMSI(int kMax, PRECISION epsilon)
{
  PartialLUDecomp(0.9);
  
  //begin iterations
  int t=0;
  for(; t<kMax; t++)
    {
      CalcResidual();
      StepMSI();
      PRECISION delta=del.Max();
      
      for(int i=0; i<=M; i++)
	{
	  for(int j=0; j<=N; j++)
	    {
	      sol[i][j]+=del[i][j];
	    }
	}

      if(fabs(delta)<epsilon) break;
    }

  return t;
}

int TNonaDiagMatrix::SolveMSI(int kMax, PRECISION epsilon)
{
  PartialLUDecomp(0.9);
  
  int t=0;
  for(; t<kMax; t++)
    {
      CalcResidual();
      StepMSI();
      PRECISION delta=del.Max();
      
      for(int i=0; i<=M; i++)
	{
	  for(int j=0; j<=N; j++)
	    {
	      sol[i][j]+=del[i][j];
	    }
	}

      if(fabs(delta)<epsilon) break;
    }

  return t;
}

int T13DiagMatrix::SolveMSI(int kMax, PRECISION epsilon)
{
  PartialLUDecomp(0.5);

  int t=0;
  for(; t<kMax; t++)
    {
      CalcResidual();
      StepMSI();
      PRECISION tError=del.Max();

      for(int i=0; i<=M; i++)
	{
	  for(int j=0; j<=N; j++)
	    {
	      sol[i][j]+=del[i][j];
	    }
	}
      
      //      cout<<tError<<endl;
      if(fabs(tError)<epsilon) break;
    }

  return t;
}

int TPentaDiagMatrix::StepMSI(void)
  /*solves the equation res=LU(del) for del by backward-forward substitution
   */
{
  int jl=0;
  int ju=N;
  int il=0;
  int iu=M;

  for(int i=il; i<=iu; i++)
    {
      register const int im=i-1;
      for(int j=jl; j<=ju; j++)
	{
	  del[i][j]=(res[i][j]-B[i][j]*del[im][j]
		     -C[i][j]*del[im][j+1]-D[i][j]*del[i][j-1])*E[i][j];
	}
    }
  
  for(int i=iu; i>=il; i--)
    {
      register const int ip=i+1;
      for(int j=ju; j>=jl; j--)
	{
	  del[i][j]-=F[i][j]*del[i][j+1]+G[i][j]*del[ip][j-1]
	    +H[i][j]*del[ip][j];
	}
    }
  
  return 1;
}

int TNonaDiagMatrix::StepMSI(void)
  /*solves the equation res=LU(del) for del by backward-forward substitution
   */
{
  int jl=0;
  int ju=N;
  int il=0;
  int iu=M;
  
  for(int i=il; i<=iu; i++)
    {
      register const int im=i-1;
      for(int j=jl; j<=ju; j++)
	{
	  del[i][j]=(res[i][j]-A[i][j]*del[im][j-1]-B[i][j]*del[im][j]
		     -C[i][j]*del[im][j+1]-D[i][j]*del[i][j-1])*E[i][j];
	}
    }
  
  for(int i=iu; i>=il; i--)
    {
      register const int ip=i+1;
      for(int j=ju; j>=jl; j--)
	{
	  del[i][j]-=F[i][j]*del[i][j+1]+G[i][j]*del[ip][j-1]
	    +H[i][j]*del[ip][j]+U[i][j]*del[ip][j+1];
	}
    }
  
  return 1;
}

int T13DiagMatrix::StepMSI(void)
  /*solves the equation res=LU(del) for del by backward-forward substitution
   */
{
  for(int i=0; i<=M; i++)
    {
      int im=i-1;
      for(int j=0; j<=N; j++)
	{
	  del[i][j]=(res[i][j]
		     -P[i][j]*del[i-2][j]
		     -A[i][j]*del[im][j-1]
		     -B[i][j]*del[im][j]
		     -C[i][j]*del[im][j+1]
		     -Q[i][j]*del[i][j-2]
		     -D[i][j]*del[i][j-1])
	    *E[i][j];
	}
    }
  
  for(int i=M; i>=0; i--)
    {
      int ip=i+1;
      for(int j=N; j>=0; j--)
	{
	  del[i][j]-=F[i][j]*del[i][j+1]
	    +R[i][j]*del[i][j+2]
	    +G[i][j]*del[ip][j-1]
	    +H[i][j]*del[ip][j]
	    +U[i][j]*del[ip][j+1]
	    +S[i][j]*del[i+2][j];
	}
    }
  
  return 1;
}

int TPentaDiagMatrix::StepGS(void)
{
  int jl=0;
  int ju=N;
  int il=0;
  int iu=M;
  
  for(int i=iu; i>=il; i--)
    {
      register const int ip=i+1;
      for(int j=ju; j>=jl; j--)
	{
	  del[i][j]=(res[i][j]
		     -a8[i][j]*del[i][j+1]
		     -a6[i][j]*del[ip][j])
	    /a5[i][j];
	}
    }
  
  return 1;
}

int TNonaDiagMatrix::StepGS(void)
{
  int jl=0;
  int ju=N;
  int il=0;
  int iu=M;
  
  for(int i=iu; i>=il; i--)
    {
      register const int ip=i+1;
      for(int j=ju; j>=jl; j--)
	{
	  del[i][j]=(res[i][j]
		     -a8[i][j]*del[i][j+1]
		     -a3[i][j]*del[ip][j-1]
		     -a6[i][j]*del[ip][j]
		     -a9[i][j]*del[ip][j+1])
	    /a5[i][j];
	}
    }
  
  return 1;
}

int T13DiagMatrix::StepGS(void)
{
  int jl=0;
  int ju=N;
  int il=0;
  int iu=M;
  
  for(int i=iu; i>=il; i--)
    {
      register const int ip=i+1;
      for(int j=ju; j>=jl; j--)
	{
	  del[i][j]=(res[i][j]
		     -a8[i][j]*del[i][j+1]
		     -ann[i][j]*del[i][j+2]
		     -a3[i][j]*del[ip][j-1]
		     -a6[i][j]*del[ip][j]
		     -a9[i][j]*del[ip][j+1]
		     -aee[i][j]*del[i+2][j])
	    /a5[i][j];
	}
    }
  
  return 1;
}

PRECISION TPentaDiagMatrix::RelaxRedBlack(int kMax, PRECISION epsilon)
{
  PRECISION resmax=1E9;
  for(int t=0; t<kMax; t++)
    {
      resmax=0;
      //use red-black relaxation
      for(int i=0; i<=M; i++)
      	{
	  for(int j=(i%2 ? 0 : 1); j<=N; j+=2)
	    {
	      PRECISION s_new=(rhs[i][j]
			       -a2[i][j]*sol[i][j-1]
			       -a4[i][j]*sol[i-1][j]
			       -a6[i][j]*sol[i+1][j]
			       -a8[i][j]*sol[i][j+1])
		/a5[i][j];
	      
	      if(fabs(s_new-sol[i][j])>resmax) resmax=fabs(s_new-sol[i][j]);
	      sol[i][j]=s_new;
            }
	}
      
      for(int i=0; i<=M; i++)
      	{
	  for(int j=(i%2 ? 1 : 0); j<=N; j+=2)
	    {
	      PRECISION s_new=(rhs[i][j]
			       -a2[i][j]*sol[i][j-1]
			       -a4[i][j]*sol[i-1][j]
			       -a6[i][j]*sol[i+1][j]
			       -a8[i][j]*sol[i][j+1])
		/a5[i][j];
	      
	      if(fabs(s_new-sol[i][j])>resmax) resmax=fabs(s_new-sol[i][j]);
	      sol[i][j]=s_new;
            }
	}

      if(resmax<epsilon) break;
    }
  
  return resmax;
}

PRECISION TNonaDiagMatrix::RelaxRedBlack(int kMax, PRECISION epsilon)
{
  PRECISION resmax=1E9;
  for(int t=0; t<kMax; t++)
    {
      resmax=0;
      //use red-black relaxation
      for(int i=0; i<=M; i++)
      	{
	  for(int j=(i%2 ? 0 : 1); j<=N; j+=2)
	    {
	      PRECISION s_new=(rhs[i][j]
			       -a1[i][j]*sol[i-1][j-1]
			       -a2[i][j]*sol[i][j-1]
			       -a3[i][j]*sol[i+1][j-1]
			       -a4[i][j]*sol[i-1][j]
			       -a6[i][j]*sol[i+1][j]
			       -a7[i][j]*sol[i-1][j+1]
			       -a8[i][j]*sol[i][j+1]
			       -a9[i][j]*sol[i+1][j+1])
		/a5[i][j];
	      
	      if(fabs(s_new-sol[i][j])>resmax) resmax=fabs(s_new-sol[i][j]);
	      sol[i][j]=s_new;
            }
	}
      
      for(int i=0; i<=M; i++)
      	{
	  for(int j=(i%2 ? 1 : 0); j<=N; j+=2)
	    {
	      PRECISION s_new=(rhs[i][j]
			       -a1[i][j]*sol[i-1][j-1]
			       -a2[i][j]*sol[i][j-1]
			       -a3[i][j]*sol[i+1][j-1]
			       -a4[i][j]*sol[i-1][j]
			       -a6[i][j]*sol[i+1][j]
			       -a7[i][j]*sol[i-1][j+1]
			       -a8[i][j]*sol[i][j+1]
			       -a9[i][j]*sol[i+1][j+1])
		/a5[i][j];
	      
	      if(fabs(s_new-sol[i][j])>resmax) resmax=fabs(s_new-sol[i][j]);
	      sol[i][j]=s_new;
            }
	}

      if(resmax<epsilon) break;
    }
  
  return resmax;
}

PRECISION T13DiagMatrix::RelaxRedBlack(int kMax, PRECISION epsilon)
{
  PRECISION resmax=1E9;
  for(int t=0; t<kMax; t++)
    {
      resmax=0;
      //use red-black relaxation
      for(int i=0; i<=M; i++)
      	{
	  for(int j=(i%2 ? 0 : 1); j<=N; j+=2)
	    {
	      PRECISION s_new=(rhs[i][j]
			       -ass[i][j]*sol[i][j-2]
			       -a1[i][j]*sol[i-1][j-1]
			       -a2[i][j]*sol[i][j-1]
			       -a3[i][j]*sol[i+1][j-1]
			       -aww[i][j]*sol[i-2][j]
			       -a4[i][j]*sol[i-1][j]
			       -a6[i][j]*sol[i+1][j]
			       -aee[i][j]*sol[i+2][j]
			       -a7[i][j]*sol[i-1][j+1]
			       -a8[i][j]*sol[i][j+1]
			       -a9[i][j]*sol[i+1][j+1]
			       -ann[i][j]*sol[i][j+2])
		/a5[i][j];
	      
	      if(fabs(s_new-sol[i][j])>resmax) resmax=fabs(s_new-sol[i][j]);
	      sol[i][j]=s_new;
            }
	}
      
      for(int i=0; i<=M; i++)
      	{
	  for(int j=(i%2 ? 1 : 0); j<=N; j+=2)
	    {
	      PRECISION s_new=(rhs[i][j]
			       -ass[i][j]*sol[i][j-2]
			       -a1[i][j]*sol[i-1][j-1]
			       -a2[i][j]*sol[i][j-1]
			       -a3[i][j]*sol[i+1][j-1]
			       -aww[i][j]*sol[i-2][j]
			       -a4[i][j]*sol[i-1][j]
			       -a6[i][j]*sol[i+1][j]
			       -aee[i][j]*sol[i+2][j]
			       -a7[i][j]*sol[i-1][j+1]
			       -a8[i][j]*sol[i][j+1]
			       -a9[i][j]*sol[i+1][j+1]
			       -ann[i][j]*sol[i][j+2])
		/a5[i][j];
	      
	      if(fabs(s_new-sol[i][j])>resmax) resmax=fabs(s_new-sol[i][j]);
	      sol[i][j]=s_new;
            }
	}

      if(resmax<epsilon) break;
    }
  
  return resmax;
}

PRECISION TNonaDiagMatrix::RelaxRedBlackOGrid(int kMax, PRECISION epsilon)
{
  PRECISION resmax=1E9;
  for(int t=0; t<kMax; t++)
    {
      resmax=0;
      //use red-black relaxation
      for(int i=0; i<=M; i++)
      	{
	  for(int j=(i%2 ? 0 : 1); j<N; j+=2)
	    {
	      int jm=j ? j-1 : N-1;
	      int jp=j+1;//j is always less than N ... =j<N ? j+1 : 1;
	      
	      PRECISION s_new=(rhs[i][j]
			       -a1[i][j]*sol[i-1][jm]
			       -a2[i][j]*sol[i][jm]
			       -a3[i][j]*sol[i+1][jm]
			       -a4[i][j]*sol[i-1][j]
			       -a6[i][j]*sol[i+1][j]
			       -a7[i][j]*sol[i-1][jp]
			       -a8[i][j]*sol[i][jp]
			       -a9[i][j]*sol[i+1][jp])
		/a5[i][j];
	      
	      if(fabs(s_new-sol[i][j])>resmax) resmax=fabs(s_new-sol[i][j]);
	      sol[i][j]=s_new;
            }
	  
	  if((i+N)%2) sol[i][N]=sol[i][0];
	}
      
      for(int i=0; i<=M; i++)
      	{
	  for(int j=(i%2 ? 1 : 0); j<N; j+=2)
	    {
	      int jm=j ? j-1 : N-1;
	      int jp=j+1;//j is always less than N ... =j<N ? j+1 : 1;
	      
	      PRECISION s_new=(rhs[i][j]
			       -a1[i][j]*sol[i-1][jm]
			       -a2[i][j]*sol[i][jm]
			       -a3[i][j]*sol[i+1][jm]
			       -a4[i][j]*sol[i-1][j]
			       -a6[i][j]*sol[i+1][j]
			       -a7[i][j]*sol[i-1][jp]
			       -a8[i][j]*sol[i][jp]
			       -a9[i][j]*sol[i+1][jp])
		/a5[i][j];
	      
	      if(fabs(s_new-sol[i][j])>resmax) resmax=fabs(s_new-sol[i][j]);
	      sol[i][j]=s_new;
            }

	  if(!((i+N)%2)) sol[i][N]=sol[i][0];
	}

      if(resmax<epsilon) break;
    }
  
  return resmax;
}

PRECISION T13DiagMatrix::RelaxRedBlackOGrid(int kMax, PRECISION epsilon)
{
  PRECISION resmax=1E9;
  for(int t=0; t<kMax; t++)
    {
      resmax=0;
      //use red-black relaxation
      for(int i=0; i<=M; i++)
      	{
	  for(int j=(i%2 ? 0 : 1); j<N; j+=2)
	    {
	      int jm=j ? j-1 : N-1;
	      int jmm=jm ? jm-1 : N-1;
	      
	      int jp=j+1;//j is always less than N ... =j<N ? j+1 : 1;
	      int jpp=jp<N ? jp+1 : 1;
	      
	      PRECISION s_new=(rhs[i][j]
			       -ass[i][j]*sol[i][jmm]
			       -a1[i][j]*sol[i-1][jm]
			       -a2[i][j]*sol[i][jm]
			       -a3[i][j]*sol[i+1][jm]
			       -aww[i][j]*sol[i-2][j]
			       -a4[i][j]*sol[i-1][j]
			       -a6[i][j]*sol[i+1][j]
			       -aee[i][j]*sol[i+2][j]
			       -a7[i][j]*sol[i-1][jp]
			       -a8[i][j]*sol[i][jp]
			       -a9[i][j]*sol[i+1][jp]
			       -ann[i][j]*sol[i][jpp])
		/a5[i][j];
	      
	      if(fabs(s_new-sol[i][j])>resmax) resmax=fabs(s_new-sol[i][j]);
	      sol[i][j]=s_new;
            }
	  
	  if((i+N)%2) sol[i][N]=sol[i][0];
	}
      
      for(int i=0; i<=M; i++)
      	{
	  for(int j=(i%2 ? 1 : 0); j<N; j+=2)
	    {
	      int jm=j ? j-1 : N-1;
	      int jmm=jm ? jm-1 : N-1;
	      
	      int jp=j+1;//j is always less than N ... =j<N ? j+1 : 1;
	      int jpp=jp<N ? jp+1 : 1;
	      
	      PRECISION s_new=(rhs[i][j]
			       -ass[i][j]*sol[i][jmm]
			       -a1[i][j]*sol[i-1][jm]
			       -a2[i][j]*sol[i][jm]
			       -a3[i][j]*sol[i+1][jm]
			       -aww[i][j]*sol[i-2][j]
			       -a4[i][j]*sol[i-1][j]
			       -a6[i][j]*sol[i+1][j]
			       -aee[i][j]*sol[i+2][j]
			       -a7[i][j]*sol[i-1][jp]
			       -a8[i][j]*sol[i][jp]
			       -a9[i][j]*sol[i+1][jp]
			       -ann[i][j]*sol[i][jpp])
		/a5[i][j];
	      
	      if(fabs(s_new-sol[i][j])>resmax) resmax=fabs(s_new-sol[i][j]);
	      sol[i][j]=s_new;
            }

	  if(!((i+N)%2)) sol[i][N]=sol[i][0];
	}

      if(resmax<epsilon) break;
    }
  
  return resmax;
}

int TPentaDiagMatrix::CalcResidual(void)
  /*Calculates the residual(?) of Ax=b. That is, it finds
    res=b-Ax' where x' is the latest estimate for x
  */
{ 
  for(int i=0; i<=M; i++)
    {
      for(int j=0; j<=N; j++)
	{
	  res[i][j]=rhs[i][j]
	    -a2[i][j]*sol[i][j-1]
	    -a4[i][j]*sol[i-1][j]
	    -a5[i][j]*sol[i][j]
	    -a6[i][j]*sol[i+1][j]
	    -a8[i][j]*sol[i][j+1];
	}
    }
  
  return 1;
}	

int TNonaDiagMatrix::CalcResidual(void)
  /*Calculates the residual(?) of Ax=b. That is, it finds
    res=b-Ax' where x' is the latest estimate for x
  */
{ 
  for(int i=0; i<=M; i++)
    {
      for(int j=0; j<=N; j++)
	{
	  res[i][j]=rhs[i][j]
	    -a1[i][j]*sol[i-1][j-1]
	    -a2[i][j]*sol[i][j-1]
	    -a3[i][j]*sol[i+1][j-1]
	    -a4[i][j]*sol[i-1][j]
	    -a5[i][j]*sol[i][j]
	    -a6[i][j]*sol[i+1][j]
	    -a7[i][j]*sol[i-1][j+1]
	    -a8[i][j]*sol[i][j+1]
	    -a9[i][j]*sol[i+1][j+1];
	}
    }
  
  return 1;
}	

int T13DiagMatrix::CalcResidual(void)
  /*Calculates the residual(?) of Ax=b. That is, it finds
    res=b-Ax' where x' is the latest estimate for x
  */
{ 
  for(int i=0; i<=M; i++)
    {
      for(int j=0; j<=N; j++)
	{
	  res[i][j]=rhs[i][j]
	    -ass[i][j]*sol[i][j-2]
	    -a1[i][j]*sol[i-1][j-1]
	    -a2[i][j]*sol[i][j-1]
	    -a3[i][j]*sol[i+1][j-1]
	    -aww[i][j]*sol[i-2][j]
	    -a4[i][j]*sol[i-1][j]
	    -a5[i][j]*sol[i][j]
	    -a6[i][j]*sol[i+1][j]
	    -aee[i][j]*sol[i+2][j]
	    -a7[i][j]*sol[i-1][j+1]
	    -a8[i][j]*sol[i][j+1]
	    -a9[i][j]*sol[i+1][j+1]
	    -ann[i][j]*sol[i][j+2];
	}
    }
  
  return 1;
}	

int TPentaDiagMatrix::Output(char* zFile)
{
  ofstream outfile(zFile);
  outfile<<"variables=\"a2\" \"a4\" \"a5\" \"a6\" \"a8\" \"rhs\" \"sol\"\"B\" \"C\" \"D\" \"E\" \"F\" \"G\" \"H\""<<endl;
  for(int i=0; i<=M; i++)
    {
      for(int j=0; j<=N; j++)
	{
	  outfile
	    <<i<<" "
	    <<j<<" "
	    <<a2[i][j]<<" "
	    <<a4[i][j]<<" "
	    <<a5[i][j]<<" "
	    <<a6[i][j]<<" "
	    <<a8[i][j]<<" "
	    <<rhs[i][j]<<" "
	    <<sol[i][j]<<" "
	    <<res[i][j]<<" "

	    <<B[i][j]<<" "
	    <<C[i][j]<<" "
	    <<D[i][j]<<" "
	    <<E[i][j]<<" "
	    <<F[i][j]<<" "
	    <<G[i][j]<<" "
	    <<H[i][j]<<" "

	    <<endl;
	}
    }

  return 1;
}

int TNonaDiagMatrix::Output(char* zFile)
{
  ofstream outfile(zFile);
  outfile<<"variables=\"a1\" \"a2\" \"a3\" \"a4\" \"a5\" \"a6\" \"a7\" \"a8\" \"a9\" \"rhs\" \"s\" \"res\" \"A\" \"B\" \"C\" \"D\" \"E\" \"F\" \"G\" \"H\" \"U\""<<endl;
  for(int i=0; i<=M; i++)
    {
      for(int j=0; j<=N; j++)
	{
	  outfile
	    <<a1[i][j]<<" "
	    <<a2[i][j]<<" "
	    <<a3[i][j]<<" "
	    <<a4[i][j]<<" "
	    <<a5[i][j]<<" "
	    <<a6[i][j]<<" "
	    <<a7[i][j]<<" "
	    <<a8[i][j]<<" "
	    <<a9[i][j]<<" "
	    <<rhs[i][j]<<" "
	    <<sol[i][j]<<" "
	    <<res[i][j]<<" "

	    <<A[i][j]<<" "
	    <<B[i][j]<<" "
	    <<C[i][j]<<" "
	    <<D[i][j]<<" "
	    <<E[i][j]<<" "
	    <<F[i][j]<<" "
	    <<G[i][j]<<" "
	    <<H[i][j]<<" "
	    <<U[i][j]

	    <<endl;
	}
    }

  return 1;
}

int T13DiagMatrix::Output(char* zFile)
{
  ofstream outfile(zFile);
  outfile<<"variables=\"x\" \"y\" \"a1\" \"a2\" \"a3\" \"a4\" \"a5\" \"a6\" \"a7\" \"a8\" \"a9\" \"rhs\" \"s\" \"A\" \"B\" \"C\" \"D\" \"E\" \"F\" \"G\" \"H\" \"U\""<<endl;
  outfile<<"zone i="<<N-1<<" j="<<M-1<<endl;
  for(int i=1; i<M; i++)
    {
      for(int j=1; j<N; j++)
	{
	  outfile<<i/(float)(M-2)<<" "<<j/(float)(N-2)<<" "
		 <<a1[i][j]<<" "
		 <<a2[i][j]<<" "
		 <<a3[i][j]<<" "
		 <<a4[i][j]<<" "
		 <<a5[i][j]<<" "
		 <<a6[i][j]<<" "
		 <<a7[i][j]<<" "
		 <<a8[i][j]<<" "
		 <<a9[i][j]<<" "

		 <<rhs[i][j]<<" "
		 <<sol[i][j]<<" "

		 <<A[i][j]<<" "
		 <<B[i][j]<<" "
		 <<C[i][j]<<" "
		 <<D[i][j]<<" "
		 <<E[i][j]<<" "
		 <<F[i][j]<<" "
		 <<G[i][j]<<" "
		 <<H[i][j]<<" "
		 <<U[i][j]

		 <<endl;
	}
    }

  return 1;
}

int TPentaDiagMatrix::PartialLUDecompOGrid(float alpha)
{
  //temporarily zero coefficients on periodic boundary
  TVector<PRECISION> temp2(M+1);
  TVector<PRECISION> temp8(M+1);

  for(int i=0; i<=M; i++)
    {
      temp2[i]=a2[i][0];
      a2[i][0]=0;
      temp8[i]=a8[i][N];
      a8[i][N]=0;
    }
  
  //l-u decomposition
  PartialLUDecomp(alpha);

  //replace the coefficients
  for(int i=0; i<=M; i++)
    {
      a2[i][0]=temp2[i];
      a8[i][N]=temp8[i];
    }
  
  return 1;
}

int TPentaDiagMatrix::SolveMSIOGrid(int kMax, PRECISION epsilon)
{
  PartialLUDecompOGrid(0.9);

  //begin iterations
  int t=0;
  for(; t<kMax; t++)
    {
      CalcResidualOGrid();
      
      StepMSI();
      PRECISION delta=del.Max();
      
      for(int i=0; i<=M; i++)
	{
	  for(int j=0; j<=N; j++)
	    {
	      sol[i][j]+=del[i][j];
	    }
	}

      if(fabs(delta)<epsilon) break;
    }

  return t;
}

int TPentaDiagMatrix::CalcResidualOGrid(void)
  /*Calculates the residual of Ax=b. That is, it finds
    res=b-Ax' where x' is the latest estimate for x
    
    assumes an o-grid that is periodic in j
  */
{ 
  for(int i=0; i<=M; i++)
    {
      for(int j=0; j<=N; j++)
	{
	  int jm=j ? j-1 : N-1;
	  int jp=j<N ? j+1 : 1;

	  res[i][j]=rhs[i][j]
	    -a2[i][j]*sol[i][jm]
	    -a4[i][j]*sol[i-1][j]
	    -a5[i][j]*sol[i][j]
	    -a6[i][j]*sol[i+1][j]
	    -a8[i][j]*sol[i][jp];
	}
    }
  
  return 1;
}	

int TNonaDiagMatrix::PartialLUDecompOGrid(float alpha)
{
  //temporarily zero coefficients on periodic boundary
  TVector<PRECISION> temp1(M+1);
  TVector<PRECISION> temp2(M+1);
  TVector<PRECISION> temp3(M+1);
  TVector<PRECISION> temp7(M+1);
  TVector<PRECISION> temp8(M+1);
  TVector<PRECISION> temp9(M+1);

  for(int i=0; i<=M; i++)
    {
      temp1[i]=a1[i][0];
      a1[i][0]=0;
      temp2[i]=a2[i][0];
      a2[i][0]=0;
      temp3[i]=a3[i][0];
      a3[i][0]=0;
      temp7[i]=a7[i][N];
      a7[i][N]=0;
      temp8[i]=a8[i][N];
      a8[i][N]=0;
      temp9[i]=a9[i][N];
      a9[i][N]=0;
    }
  
  //l-u decomposition
  PartialLUDecomp(alpha);

  //replace the coefficients
  for(int i=0; i<=M; i++)
    {
      a1[i][0]=temp1[i];
      a2[i][0]=temp2[i];
      a3[i][0]=temp3[i];
      a7[i][N]=temp7[i];
      a8[i][N]=temp8[i];
      a9[i][N]=temp9[i];
    }
  
  return 1;
}

int TNonaDiagMatrix::SolveMSIOGrid(int kMax, PRECISION epsilon)
{
  PartialLUDecompOGrid(0.9);

  //begin iterations
  int t=0;
  for(; t<kMax; t++)
    {
      CalcResidualOGrid();
      
      StepMSI();
      PRECISION delta=del.Max();
      
      for(int i=0; i<=M; i++)
	{
	  for(int j=0; j<=N; j++)
	    {
	      sol[i][j]+=del[i][j];
	    }
	}

      if(fabs(delta)<epsilon) break;
    }

  return t;
}

int TNonaDiagMatrix::CalcResidualOGrid(void)
  /*Calculates the residual of Ax=b. That is, it finds
    res=b-Ax' where x' is the latest estimate for x
    
    assumes an o-grid that is periodic in j
  */
{ 
  for(int i=0; i<=M; i++)
    {
      for(int j=0; j<=N; j++)
	{
	  int jm=j ? j-1 : N-1;
	  int jp=j<N ? j+1 : 1;

	  res[i][j]=rhs[i][j]
	    -a1[i][j]*sol[i-1][jm]
	    -a2[i][j]*sol[i][jm]
	    -a3[i][j]*sol[i+1][jm]
	    -a4[i][j]*sol[i-1][j]
	    -a5[i][j]*sol[i][j]
	    -a6[i][j]*sol[i+1][j]
	    -a7[i][j]*sol[i-1][jp]
	    -a8[i][j]*sol[i][jp]
	    -a9[i][j]*sol[i+1][jp];
	}
    }
  
  return 1;
}	

int T13DiagMatrix::PartialLUDecompOGrid(float alpha)
{
  //temporarily zero coefficients on periodic boundary
  TVector<PRECISION> tempss(M+1);
  TVector<PRECISION> tempss1(M+1);
  TVector<PRECISION> temp1(M+1);
  TVector<PRECISION> temp2(M+1);
  TVector<PRECISION> temp3(M+1);
  TVector<PRECISION> temp7(M+1);
  TVector<PRECISION> temp8(M+1);
  TVector<PRECISION> temp9(M+1);
  TVector<PRECISION> tempnn(M+1);
  TVector<PRECISION> tempnnN1(M+1);

  for(int i=0; i<=M; i++)
    {
      tempss1[i]=ass[i][1];
      ass[i][1]=0;
      tempss[i]=ass[i][0];
      ass[i][0]=0;
      temp1[i]=a1[i][0];
      a1[i][0]=0;
      temp2[i]=a2[i][0];
      a2[i][0]=0;
      temp3[i]=a3[i][0];
      a3[i][0]=0;
      temp7[i]=a7[i][N];
      a7[i][N]=0;
      temp8[i]=a8[i][N];
      a8[i][N]=0;
      temp9[i]=a9[i][N];
      a9[i][N]=0;
      tempnn[i]=ann[i][N];
      ann[i][N]=0;
      tempnnN1[i]=ann[i][N-1];
      ann[i][N-1]=0;
    }
  
  //l-u decomposition
  PartialLUDecomp(alpha);

  //replace the coefficients
  for(int i=0; i<=M; i++)
    {
      ass[i][1]=tempss1[i];
      ass[i][0]=tempss[i];
      a1[i][0]=temp1[i];
      a2[i][0]=temp2[i];
      a3[i][0]=temp3[i];
      a7[i][N]=temp7[i];
      a8[i][N]=temp8[i];
      a9[i][N]=temp9[i];
      ann[i][N]=tempnn[i];
      ann[i][N-1]=tempnnN1[i];
    }
  
  return 1;
}

PRECISION T13DiagMatrix::SolveMSIOGrid(int kMax, PRECISION epsilon)
{
  //begin iterations
  int t=0;
  PRECISION delta=0;
  
  for(; t<kMax; t++)
    {
      CalcResidualOGrid();
      
      StepMSI();
      if(epsilon>0) delta=del.Max();
      
      for(int i=0; i<=M; i++)
	{
	  for(int j=0; j<N; j++)
	    {
	      sol[i][j]+=del[i][j];
	    }
	  sol[i][N]=sol[i][0];
	}

      if(fabs(delta)<epsilon) break;
    }

  return delta;
}

int T13DiagMatrix::CalcResidualOGrid(void)
  /*Calculates the residual of Ax=b. That is, it finds
    res=b-Ax' where x' is the latest estimate for x
    
    assumes an o-grid that is periodic in j
  */
{ 
  for(int i=0; i<=M; i++)
    {
      for(int j=0; j<N; j++)
	{
	  int jm=j ? j-1 : N-1;
	  int jmm=jm ? jm-1 : N-1;

	  int jp=j+1;//j is always less than N
	  int jpp=jp<N ? jp+1 : 1;

	  res[i][j]=rhs[i][j]
	    -ass[i][j]*sol[i][jmm]
  	    //-a1[i][j]*sol[i-1][jm]
	    -a2[i][j]*sol[i][jm]
  	    //-a3[i][j]*sol[i+1][jm]
	    -aww[i][j]*sol[i-2][j]
	    -a4[i][j]*sol[i-1][j]
	    -a5[i][j]*sol[i][j]
	    -a6[i][j]*sol[i+1][j]
	    -aee[i][j]*sol[i+2][j]
  	    //-a7[i][j]*sol[i-1][jp]
	    -a8[i][j]*sol[i][jp]
  	    //-a9[i][j]*sol[i+1][jp]
	    -ann[i][j]*sol[i][jpp];
	}

      res[i][N]=res[i][0];
    }
  
  return 1;
}	

int TPentaDiagMatrix::SetupMSI(void)
{
  B=matrix(M+5, N+5, -2, -2);
  C=matrix(M+5, N+5, -2, -2);
  D=matrix(M+5, N+5, -2, -2);
  E=matrix(M+5, N+5, -2, -2);
  F=matrix(M+5, N+5, -2, -2);
  G=matrix(M+5, N+5, -2, -2);
  H=matrix(M+5, N+5, -2, -2);
  
  return 1;
}

int TNonaDiagMatrix::SetupMSI(void)
{
  TPentaDiagMatrix::SetupMSI();
  A=matrix(M+5, N+5, -2, -2);
  U=matrix(M+5, N+5, -2, -2);
  
  return 1;
}

int T13DiagMatrix::SetupMSI(void)
{
  TNonaDiagMatrix::SetupMSI();
  
  P=matrix(M+5, N+5, -2, -2);
  Q=matrix(M+5, N+5, -2, -2);
  R=matrix(M+5, N+5, -2, -2);
  S=matrix(M+5, N+5, -2, -2);
  
  return 1;
}

TBanded3D::TBanded3D(int I, int J, int K) : x(I,J,K),
					    a0(I,J,K), aW(I,J,K), aE(I,J,K), aS(I,J,K),
					    aN(I,J,K), aU(I,J,K), aD(I,J,K), rhs(I,J,K)//,
  //R(I,J,K)//, f(I,J,K), R0(I,J,K), def(I,J,K)
{
  iMax=I-1; jMax=J-1; kMax=K-1;

  for(int i=0; i<I; i++)
    for(int j=0; j<J; j++)
      for(int k=0; k<K; k++)
	{
	  rhs[i][j][k]=dvector(4);
	  x[i][j][k]=dvector(4);

	  a0[i][j][k]=dmatrix(4,4);
	  aW[i][j][k]=dmatrix(4,4);
	  aE[i][j][k]=dmatrix(4,4);
	  aS[i][j][k]=dmatrix(4,4);
	  aN[i][j][k]=dmatrix(4,4);
	  aD[i][j][k]=dmatrix(4,4);
	  aU[i][j][k]=dmatrix(4,4);

	  //R[i][j][k]=dvector(4);
	  /*R0[i][j][k]=dvector(4);
	    f[i][j][k]=dvector(4);
	    def[i][j][k]=dvector(4);
	  */
	}
}

int TBanded3D::RelaxGauss( int nu )
{
  T5DMatrix a0inv(iMax+1, jMax+1, kMax+1);

      for(int i=0; i<=iMax; i++)
	{
	  for(int j=0; j<=jMax; j++)
	    {
	      for(int k=0; k<=kMax; k++)
		{
		  a0inv[i][j][k]=a0[i][j][k].FindInverse();
		}
	    }
	}

      //T4DMatrix xnew=x;
  for(int pass=0; pass<nu; pass++)
    {
      for(int i=0; i<=iMax; i++)
	{
	  for(int j=0; j<=jMax; j++)
	    {
	      for(int k=0; k<=kMax; k++)
		{
		  dvector xrhs=rhs[i][j][k];

		  //really slow with all the if's, but this is a slow method anyway
		  if(i>0)  xrhs-=aW[i][j][k]*x[i-1][j][k];
		  if(i<iMax)  xrhs-=aE[i][j][k]*x[i+1][j][k];

		  if(j>0)  xrhs-=aS[i][j][k]*x[i][j-1][k];
		  if(j<jMax)  xrhs-=aN[i][j][k]*x[i][j+1][k];

		  if(k>0)  xrhs-=aD[i][j][k]*x[i][j][k-1];
		  if(k<kMax)  xrhs-=aU[i][j][k]*x[i][j][k+1];

		  x[i][j][k]=a0inv[i][j][k]*xrhs;
		}
	    }
	}
      //x=xnew;
    }

  return 0;
}

int TBanded3D::LineRelax( )
{
  for(int j=0; j<=jMax; j++)
    {
      for(int k=0; k<=kMax; k++)
	{
	  LineRelaxI(j,k);
	}
    }

  for(int i=0; i<=iMax; i++)
    {
      for(int k=0; k<=kMax; k++)
	{
	  LineRelaxJ(i,k);
	}
    }

  if(kMax)
    {
      for(int i=0; i<=iMax; i++)
	{
	  for(int j=0; j<=jMax; j++)
	    {
	      LineRelaxK(i,j);
	    }
	}
    }

  for(int j=jMax; j>=0; j--)
    {
      for(int k=kMax; k>=0; k--)
	{
	  LineRelaxI(j,k);
	}
    }

  for(int i=iMax; i>=0; i--)
    {
      for(int k=kMax; k>=0; k--)
	{
	  LineRelaxJ(i,k);
	}
    }

  if(kMax)
    {
      for(int i=iMax; i>=0; i--)
	{
	  for(int j=jMax; j>=0; j--)
	    {
	      LineRelaxK(i,j);
	    }
	}
	
    }

  return 0;
}

int TBanded3D::LineRelaxI(int j, int k )
{
  TTriDiagMatrix<dmatrix> tDiagMat(iMax+1);
  TVector<dmatrix> tRHS(iMax+1);
  TVector<dmatrix> tX(iMax+1);
  for(int i=0; i<=iMax; i++)
    {
      tDiagMat._tA[i]=dmatrix(4,4);
      tDiagMat._tB[i]=dmatrix(4,4);
      tDiagMat._tC[i]=dmatrix(4,4);

      tRHS[i]=dmatrix(4,1);
      tX[i]=dmatrix(4,1);
    }

  int jm = j ? j-1 : jMax;
  int jp = j < jMax ? j+1 : 0;

  for(int i=0; i<=iMax; i++)
    {
      tDiagMat._tB[i]=a0[i][j][k];
      tRHS[i]=rhs[i][j][k];

      if(i>0)  tDiagMat._tA[i]=aW[i][j][k];
      if(i<iMax)  tDiagMat._tC[i]=aE[i][j][k];

      tRHS[i]-=aS[i][j][k]*x[i][jm][k];
      tRHS[i]-=aN[i][j][k]*x[i][jp][k];

      if(k>0)  tRHS[i]-=aD[i][j][k]*x[i][j][k-1];
      if(k<kMax)  tRHS[i]-=aU[i][j][k]*x[i][j][k+1];
    }			       

  tDiagMat.Solve(tX, tRHS);

  for(int i=0;i<=iMax;i++)
    {
      for(int var=0; var<4; var++)
	(x[i][j][k])[var]=(tX[i])[var][0];
    }

  return 0;
}

int TBanded3D::LineRelaxJ(int i, int k )
{
  TCyclicTriDiagMatrix<dmatrix> tDiagMat(jMax+1);
  TVector<dmatrix> tRHS(jMax+1);
  TVector<dmatrix> tX(jMax+1);

  for(int j=0; j<=jMax; j++)
    {
      tDiagMat._tA[j]=dmatrix(4,4);
      tDiagMat._tB[j]=dmatrix(4,4);
      tDiagMat._tC[j]=dmatrix(4,4);

      tRHS[j]=dmatrix(4,1);
      tX[j]=dmatrix(4,1);
    }

  for(int j=0; j<=jMax; j++)
    {
      tDiagMat._tB[j]=a0[i][j][k];
      tRHS[j]=rhs[i][j][k];

      tDiagMat._tA[j]=aS[i][j][k];
      tDiagMat._tC[j]=aN[i][j][k];

      if(i>0)  tRHS[j]-=aW[i][j][k]*x[i-1][j][k];
      if(i<iMax)  tRHS[j]-=aE[i][j][k]*x[i+1][j][k];

      if(k>0)  tRHS[j]-=aD[i][j][k]*x[i][j][k-1];
      if(k<kMax)  tRHS[j]-=aU[i][j][k]*x[i][j][k+1];
    }			       

  tDiagMat.Solve(tX, tRHS);

  for(int j=0;j<=jMax;j++)
    {
      for(int var=0; var<4; var++)
	(x[i][j][k])[var]=(tX[j])[var][0];
    }

  return 0;
}

int TBanded3D::LineRelaxK(int i, int j )
{
  TTriDiagMatrix<dmatrix> tDiagMat(kMax+1);
  TVector<dmatrix> tRHS(kMax+1);
  TVector<dmatrix> tX(kMax+1);

  for(int k=0; k<=kMax; k++)
    {
      tDiagMat._tA[k]=dmatrix(4,4);
      tDiagMat._tB[k]=dmatrix(4,4);
      tDiagMat._tC[k]=dmatrix(4,4);

      tRHS[k]=dmatrix(4,1);
      tX[k]=dmatrix(4,1);
    }

  int jm = j ? j-1 : jMax;
  int jp = j < jMax ? j+1 : 0;

  for(int k=0; k<=kMax; k++)
    {
      tDiagMat._tB[k]=a0[i][j][k];
      tRHS[k]=rhs[i][j][k];

      tDiagMat._tA[k]=aD[i][j][k];
      tDiagMat._tC[k]=aU[i][j][k];

      if(i>0)  tRHS[k]-=aW[i][j][k]*x[i-1][j][k];
      if(i<iMax)  tRHS[k]-=aE[i][j][k]*x[i+1][j][k];

      tRHS[k]-=aS[i][j][k]*x[i][jm][k];
      tRHS[k]-=aN[i][j][k]*x[i][jp][k];
    }			       

  tDiagMat.Solve(tX, tRHS);

  for(int k=0;k<=kMax;k++)
    {
      for(int var=0; var<4; var++)
	(x[i][j][k])[var]=(tX[k])[var][0];
    }

  return 0;
}

int TBanded3D::Zero(void)
{
  x.Zero();
  a0.Zero();
  aW.Zero();
  aE.Zero();
  aS.Zero();
  aN.Zero();
  aD.Zero();
  aU.Zero();
  rhs.Zero();

  return 1;
}

