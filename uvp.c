#include "uvp.h"
#include "math.h"


/**
 * Determines the value of U and G according to the formula
 *
 * @f$ F_{i,j} := u_{i,j} + \delta t \left( \frac{1}{Re} \left( \left[
    \frac{\partial^2 u}{\partial x^2} \right]_{i,j} + \left[
    \frac{\partial^2 u}{\partial y^2} \right]_{i,j} \right) - \left[
    \frac{\partial (u^2)}{\partial x} \right]_{i,j} - \left[
    \frac{\partial (uv)}{\partial y} \right]_{i,j} + g_x \right) @f$
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 *
 * @f$ G_{i,j} := v_{i,j} + \delta t \left( \frac{1}{Re} \left(
   \left[ \frac{\partial^2 v}{\partial x^2}\right]_{i,j} + \left[ \frac{\partial^2 v}{\partial
                   y^2} \right]_{i,j} \right) - \left[ \frac{\partial
                   (uv)}{\partial x} \right]_{i,j} - \left[
                 \frac{\partial (v^2)}{\partial y} \right]_{i,j} + g_y
               \right) @f$
 *
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 */
void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
)
{
	for (int j=0;j<=jmax;j++){
		F[0][j]=U[0][j];
		F[imax][j]=U[imax][j];
	}

	for(int i=0;i<=imax;i++){
		G[i][0]=V[i][0];
		G[i][jmax]=V[i][jmax];
	}


	double a, b,c,d;
	double ux2,uy2,u2x,uvy;
	double vy2,vx2,v2y,uvx;

		for(int i=1; i<imax; i++){
			for (int j=1; j<=jmax; j++){
				ux2 = (U[i+1][j]-2*U[i][j]+U[i-1][j])/(dx*dx);
				uy2 = (U[i][j+1]-2*U[i][j]+U[i][j-1])/(dy*dy);
				vy2 = (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy*dy);
				vx2 = (V[i+1][j]-2*V[i][j]+V[i-1][j])/(dx*dx);

				a = ((U[i][j]+U[i+1][j])/2);
				b = ((U[i-1][j]+U[i][j])/2);
				c = ((V[i][j]+V[i][j+1])/2);
				d = ((V[i][j-1]+V[i][j])/2);

				u2x = (a*a-b*b)/dx+(alpha/dx)*(fabs(a)*((U[i][j]-U[i+1][j])/2)-
						fabs(b)*((U[i-1][j]-U[i][j])/2));

				uvy = (1/(4*dy))*((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])-
						(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j]))+
						(alpha/(4*dy))*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-
						fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j]));

				v2y = (c*c-d*d)/dy+(alpha/dy)*(fabs(c)*((V[i][j]-V[i][j+1])/2) -
						fabs(d)*((V[i][j-1]-V[i][j])/2));

				uvx = (1/(4*dx))*((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])-
						(U[i-1][j] + U[i-1][j+1])*(V[i-1][j]+V[i][j]))+
						(alpha/(4*dx))*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-
						fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j]));

				F[i][j] = U[i][j] + dt*((1/Re)*(ux2+uy2) - u2x - uvy + GX);
				G[i][j] = V[i][j] + dt*((1/Re)*(vx2+vy2) - uvx - v2y + GY);
			}
		}
}


/**
 * This operation computes the right hand side of the pressure poisson equation.
 * The right hand side is computed according to the formula
 *
 * @f$ rs = \frac{1}{\delta t} \left( \frac{F^{(n)}_{i,j}-F^{(n)}_{i-1,j}}{\delta x} + \frac{G^{(n)}_{i,j}-G^{(n)}_{i,j-1}}{\delta y} \right)  @f$
 *
 */
void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
)
{
	for(int i=1; i<(imax); i++)
			for(int j=1; j<(jmax); j++)
				RS[i][j] = 1/dt*(((F[i][j]-F[i-1][j])/dx) + ((G[i][j]-G[i][j-1])/dy));
}


/**
 * Determines the maximal time step size. The time step size is restricted
 * accordin to the CFL theorem. So the final time step size formula is given
 * by
 *
 * @f$ {\delta t} := \tau \, \min\left( \frac{Re}{2}\left(\frac{1}{{\delta x}^2} + \frac{1}{{\delta y}^2}\right)^{-1},  \frac{{\delta x}}{|u_{max}|},\frac{{\delta y}}{|v_{max}|} \right) @f$
 *
 */
void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
)
{
	//double first_term = (Re/2)*pow((1/pow(dx, 2))+(1/pow(dy, 2)), -1);
	double first_term = (Re/2)*((dx*dx)*(dy*dy))/((dx*dx)+(dy*dy));

	double umax=0, vmax=0;
	for(int i=0; i<imax; i++)
		for(int j=0; j<jmax; j++)
			if(fabs(umax) < fabs(U[i][j]))
				umax = U[i][j];

	for(int i=0; i<imax; i++)
			for(int j=0; j<jmax; j++)
				if(fabs(vmax) < fabs(V[i][j]))
					vmax = V[i][j];


	double second_term = dx / fabs(umax);
	double third_term = dy / fabs(vmax);

	double first_min = fmin(first_term, second_term);
	double second_min = fmin(first_min, third_term);

	//double my_min = fmin(first_term, second_term, third_term);
	*dt = tau * second_min;
}


/**
 * Calculates the new velocity values according to the formula
 *
 * @f$ u_{i,j}^{(n+1)}  =  F_{i,j}^{(n)} - \frac{\delta t}{\delta x} (p_{i+1,j}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 * @f$ v_{i,j}^{(n+1)}  =  G_{i,j}^{(n)} - \frac{\delta t}{\delta y} (p_{i,j+1}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 *
 * As always the index range is
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 * @image html calculate_uv.jpg
 */
void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
)
{
	for (int i=1; i<=imax-1; i++)
		for (int j=1; j<=jmax; j++)
			U[i][j]=F[i][j]-(dt/dx)*(P[i+1][j]-P[i][j]);

	for (int i=1; i<=imax; i++)
		for (int j=1; j<=jmax-1; j++)
			V[i][j]=G[i][j]-(dt/dy)*(P[i][j+1]-P[i][j]);

}
