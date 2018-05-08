#include "helper.h"
#include "visual.h"
#include "init.h"
#include <stdio.h>
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args){
	//ERROR("hi1");ERROR("hi1");
	const char * szFile = "cavity100.dat";
	//ERROR("hi1\n");
	double Re=0, UI=0, VI=0, PI=0, GX=0, GY=0, t_end=0, xlength=0, ylength=0,
		   dt=0, dx=0, dy=0, alpha=0, omg=0, tau=0, eps=0, dt_value=0, t=0;
	int imax=0, jmax=0, itermax=0, n=0;

	int read = read_parameters(szFile,&Re,&UI,&VI,&PI,&GX,&GY,&t_end,&xlength,&ylength,
	&dt,&dx,&dy,&imax,&jmax,&alpha,&omg,&tau,&itermax,&eps,&dt_value);
	//ERROR("hi1");

	if (read == 1)
	{
		printf("Parameters read successfully\n");
	}
	else
	{
		printf("there was an error while reading parameters\n");
	}

    double **U = matrix(0,imax+1,0,jmax+1);
    double **V = matrix(0,imax+1,0,jmax+1);
    double **P = matrix(0,imax+1,0,jmax+1);
 	double **RS = matrix(0,imax+1,0,jmax+1);
 	double **F = matrix(0,imax+1,0,jmax+1);
 	double **G = matrix(0,imax+1,0,jmax+1);

 	t=0;
	n=0;
	int n1 = 0;
	double res;
	int it;
	init_uvp(UI,VI,PI,imax,jmax,U,V,P);
	while (t<t_end)
	{
		calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,U,V);
		//printf("%f \n",dt);
		boundaryvalues(imax,jmax,U,V);
		calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G);
		calculate_rs(dt,dx,dy,imax,jmax,F,G,RS);
		it=0;
		res = 1;
		//printf("%f \n",t);
		while (it<itermax && res>eps)
		{
			sor(omg,dx,dy,imax,jmax,P,RS,&res);
			it=it+1;
		}

		calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P);
		if(t>=n1*dt_value){
			write_vtkFile("vis", n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
			n1 = n1 + 1;
		}

		t=t+dt;
		n=n+1;
		printf("%f \n",t);
	}

	free_matrix(P, 0, imax+1, 0, jmax+1);
	free_matrix(U, 0, imax+1, 0, jmax+1);
	free_matrix(V, 0, imax+1, 0, jmax+1);
	free_matrix(F, 0, imax+1, 0, jmax+1);
	free_matrix(G, 0, imax+1, 0, jmax+1);
	free_matrix(RS, 0, imax+1, 0, jmax+1);
//printf("%f \n",t);
	ERROR("END");
	//return -1;
}
