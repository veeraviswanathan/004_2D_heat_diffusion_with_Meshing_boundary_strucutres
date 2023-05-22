#include<stdio.h>
#include<stdlib.h>
#include<array>
#include<cmath>
#include"user_header.h"


double Flux_calc(double flux_in, double flux_out,double dx);

int main()
{
	//declare domain variables
	double end_time=300 ; // seconds
	double time_now=0.0;
	double dt=1e-5; // time step
	int n_nodes_x = 50; // no of grid points
	int n_nodes_y = 50; // no of grid points
	double length_x = 0.05;// length of domain (m)
	double length_y = 0.05;// length of domain (m)
	double delta_x = (double)length_x / n_nodes_x;
	double delta_y = (double)length_y / n_nodes_y;
	
	double Tl = 600, Tr = 273; //left and right boundary condition (K)
	double Tinit = 300; // initial temperature for domain
	double k = 2e-3; // Diffusivity (m2/s)

	// declare index variables
	int i,j,indP;
	int iter = 0;
	char fname[80];

	cellcenter** cen =nullptr;
	cen=(cellcenter**)malloc(n_nodes_x*sizeof(cellcenter*));
	for(j=0;j<n_nodes_x;j++)
	{
		cen[j]=(cellcenter*)malloc(n_nodes_y*sizeof(cellcenter));
	}	

	double rhs_coefficient=k*dt;
	FILE *fp;

	// initialisation block

	for (i = 0; i< n_nodes_x;i++)
	{
		for (j = 0; j<n_nodes_y;j++)
		{
			cen[i][j].T = Tinit;
			cen[i][j].Told = cen[i][j].T;
			cen[i][j].S = 0 * (dt);
			cen[i][j].FE=0;
			cen[i][j].FW=0;
			cen[i][j].FN=0;
			cen[i][j].FS=0;
			cen[i][j].bound_id=999;
			cen[i][j].corner_flag=0;

		}

	}

	//Mesh generator
	//***************************************************//
	for (j = 0; j< n_nodes_y;j++)
	{
		cen[0][j].bound_id=1; //LEFT
		cen[n_nodes_x-1][j].bound_id=2; //RIGHT
		for (i = 0; i<n_nodes_x;i++)
		{
			cen[i][j].xcen[0]=i*length_x/(n_nodes_x-1);
			cen[i][j].xcen[1]=j*length_y/(n_nodes_y-1);		
			cen[i][0].bound_id=4;	//BOTTOM
			cen[i][n_nodes_y-1].bound_id=3; //TOP
			if(i==0 && j==0)
			{
				cen[i][j].corner_flag=1;
			}
			if(i==n_nodes_x-1 && j==0)
			{
				cen[i][j].corner_flag=1;
			}
			if(i==0 && j==n_nodes_y-1)
			{
				cen[i][j].corner_flag=1;
			}
			if(i==n_nodes_x-1 && j==n_nodes_y-1)
			{
				cen[i][j].corner_flag=1;
			}
		}
	}

	
	//boundary id allocator
	sprintf(fname, "mesh_boundary_id.dat");	
	fp=fopen(fname,"w");
	for(j=0;j<n_nodes_y;j++)
	{
		for (i = 0;i < n_nodes_x;i++)
		{
			fprintf(fp,"%d \t",cen[i][j].bound_id);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);	

 	/****************************************************/
	//Mesh generator END
	//Gauss Siedel iterative method

	while(time_now<end_time)
	{
		printf("Iteration: %d\n",iter);
		// Assign old array
		for (i = 0;i < n_nodes_x;i++)
		{
			for (j = 0;j < n_nodes_y;j++)
			{		
				cen[i][j].Told = cen[i][j].T;
			}			
		}
	
		//Main Stencil
		for (i = 0; i < n_nodes_x; i++)
		{
			for(j=0; j<n_nodes_y;j++)
			{			
				if(cen[i][j].bound_id==999)
				{
					cen[i][j].FE=Flux_calc(cen[i][j].Told,cen[i+1][j].Told,delta_x);
					cen[i][j].FW=Flux_calc(cen[i-1][j].Told,cen[i][j].Told,delta_x);
					cen[i][j].FN=Flux_calc(cen[i][j].Told,cen[i][j+1].Told,delta_y);
					cen[i][j].FS=Flux_calc(cen[i][j-1].Told,cen[i][j].Told,delta_y);				
				}
				if(cen[i][j].bound_id==1)
				{					
					cen[i][j].FE=Flux_calc(cen[i][j].Told,cen[i+1][j].Told,delta_x);
					cen[i][j].FW=Flux_calc(Tl,cen[i][j].Told,delta_x);					
					cen[i][j].FN=Flux_calc(cen[i][j].Told,cen[i][j+1].Told,delta_y);
					cen[i][j].FS=Flux_calc(cen[i][j-1].Told,cen[i][j].Told,delta_y);				

				}
				if(cen[i][j].bound_id==2)
				{
					cen[i][j].FE=Flux_calc(cen[i][j].Told,Tr,delta_x);
					cen[i][j].FW=Flux_calc(cen[i-1][j].Told,cen[i][j].Told,delta_x);
					cen[i][j].FN=Flux_calc(cen[i][j].Told,cen[i][j+1].Told,delta_y);
					cen[i][j].FS=Flux_calc(cen[i][j-1].Told,cen[i][j].Told,delta_y);
					
				}
				if(cen[i][j].bound_id==4)
				{					
					if(cen[i][j].corner_flag==1)
					{
						if(i==0)
						{
							cen[i][j].FW=Flux_calc(Tl,cen[i][j].Told,delta_x);
						}
						else if(i==n_nodes_x-1)
						{
							cen[i][j].FE=Flux_calc(cen[i][j].Told,Tr,delta_x);
						}
					}
					else
					{				
						cen[i][j].FE=Flux_calc(cen[i][j].Told,cen[i+1][j].Told,delta_x);
						cen[i][j].FW=Flux_calc(cen[i-1][j].Told,cen[i][j].Told,delta_x);
					}
					cen[i][j].FS=0;
					cen[i][j].FN=Flux_calc(cen[i][j].Told,cen[i][j+1].Told,delta_y);											
				}
				if(cen[i][j].bound_id==3)
				{
					
					if(cen[i][j].corner_flag==1)
					{
						if(i==0)
						{
							cen[i][j].FW=Flux_calc(Tl,cen[i][j].Told,delta_x);
						}
						else if(i==n_nodes_x-1)
						{
							cen[i][j].FE=Flux_calc(cen[i][j].Told,Tr,delta_x);
						}
					}
					else
					{				
						cen[i][j].FE=Flux_calc(cen[i][j].Told,cen[i+1][j].Told,delta_x);
						cen[i][j].FW=Flux_calc(cen[i-1][j].Told,cen[i][j].Told,delta_x);
					}					
					cen[i][j].FN=0;
					cen[i][j].FS=Flux_calc(cen[i][j-1].Told,cen[i][j].Told,delta_y);		
				}
				cen[i][j].T = cen[i][j].Told + (rhs_coefficient/delta_x)*(cen[i][j].FW-cen[i][j].FE)+(rhs_coefficient/delta_y)*(cen[i][j].FS-cen[i][j].FN)+dt*cen[i][j].S;			
			}
		}
			
			
		
		//Write output file		
   		sprintf(fname, "user_heat_diffusion_transient_explicit_%f.out", time_now);			
		if(iter%1000==0 || iter==1)
		{
			fp=fopen(fname,"w");
			for(j=0;j<n_nodes_y;j++)
			{
				for (i = 0;i < n_nodes_x;i++)
				{
					fprintf(fp,"%e \t",cen[i][j].T);
				}
				fprintf(fp,"\n");
			}
			fclose(fp);	
		}	
	iter++;
	time_now=time_now+iter*dt;
	}	
}
double Flux_calc(double flux_in, double flux_out, double dx)
{
	double flux;
	flux= (flux_out-flux_in)/dx;
    return -1*flux;
}






















	
// // 	//S[int(n_nodes/2)]= 5000 * (delta_x * delta_x / k);	
	
// // 	//Gauss Siedel iterative method
// // 	while(time_now<end_time)
// // 	{
		
// // 		printf("Iteration: %d\n",iter);
// // 		// Assign old array
// // 		for (i = 0;i < n_nodes_x;i++)
// // 		{
// // 			for (j = 0;j < n_nodes_x;j++)
// // 			{		
// // 				cen[i][j].Told = cen[i][j].T;
// // 			}			
// // 		}

	
// // 		// Main Stencil
// // 		for (i = 1; i < n_nodes_x-1; i++)
// // 		{
// // 			for(j=1; j<n_nodes_y-1;j++)
// // 			{				
// // 				cen[i][j].FE=Flux_calc(cen[i][j].Told,cen[i+1][j].Told,delta_x);
// // 				cen[i][j].FW=Flux_calc(cen[i-1][j].Told,cen[i][j].Told,delta_x);
// // 				cen[i][j].FN=Flux_calc(cen[i][j].Told,cen[i][j+1].Told,delta_y);
// // 				cen[i][j].FS=Flux_calc(cen[i][j-1].Told,cen[i][j].Told,delta_y);				
// // 				cen[i][j].T = cen[i][j].Told + (rhs_coefficient/delta_x)*(cen[i][j].FW-cen[i][j].FE)+(rhs_coefficient/delta_y)*(cen[i][j].FS-cen[i][j].FN)+dt*cen[i][j].S;
// // 			}
// // 		}
		
// // 		printf("Iteration: %d\n",iter);
// // 		//Write output file
// // 		//if(fmod(time_now,2.0)<5e-3)
		
// //    		sprintf(fname, "user_heat_diffusion_transient_explicit_%f.out", time_now);	
		
// // 		if(iter%1000==0)
// // 		{
// // 			fp=fopen(fname,"w");
// // 			for(j=0;j<n_nodes_y;j++)
// // 			{
// // 				for (i = 0;i < n_nodes_x;i++)
// // 				{
// // 					fprintf(fp,"%e \t",cen[i][j].T);
// // 				}
// // 				fprintf(fp,"\n");
// // 			}
// // 			fclose(fp);	
// // 		}
		
// // 		iter++;
// // 		time_now=time_now+iter*dt;
// // 	}
	
// // 	free(cen);
// // 	return(0);
// // }



			
// // 		// Boundary conditions for the rectangular domain; LEFT WALL= Tl, RIGHT wall=TR; TOP and bottom wall , Flux=0;
// // 		// //TOP and bottom
// // 		// for(i=0;i<n_nodes_x;i++)
// // 		// {
// // 		// 	cen[i][0].FN=Flux_calc(cen[i][0].Told,cen[i][1].Told,delta_y);
// // 		// 	cen[i][0].FS=0.0;
// // 		// 	if(i==0)
// // 		// 	{
// // 		// 		cen[i][0].FW=Flux_calc(2*Tl,2*cen[i][0].Told,delta_x);
// // 		// 		// cen[i][j].FW=0.0;
// // 		// 	}
// // 		// 	else
// // 		// 	{
// // 		// 		cen[i][0].FW=Flux_calc(cen[i-1][j].Told,cen[i][0].Told,delta_x);
// // 		// 	}
// // 		// 	if(i==n_nodes_x-1)
// // 		// 	{
// // 		// 		// cen[i][j].FE=Flux_calc(2*Told[i][0],2*Tr,delta_x);
// // 		// 		cen[i][0].FE=0.0;
// // 		// 	}
// // 		// 	else
// // 		// 	{
// // 		// 		cen[i][0].FE=Flux_calc(cen[i][0].Told,cen[i+1][0].Told,delta_x);
// // 		// 	}			
			
// // 		// 	// BOTTOM WALL TEMPERATURE			
// // 		// 	T[i][0] = Told[i][0] + (rhs_coefficient/delta_x)*(cen[i][j].FW-cen[i][j].FE)+(rhs_coefficient/delta_y)*(FS-FN)+dt*S[i][0];
			
// // 		// 	// #################################//
// // 		// 	FN=0.0;
// // 		// 	FS=Flux_calc(Told[i][n_nodes_y-2],Told[i][n_nodes_y-1],delta_y);				;
// // 		// 	if(i==0)
// // 		// 	{
// // 		// 		cen[i][j].FW=Flux_calc(2*Tl,2*Told[i][n_nodes_y-1],delta_x);
// // 		// 	}
// // 		// 	else
// // 		// 	{
// // 		// 		cen[i][j].FW=Flux_calc(Told[i-1][n_nodes_y-1],Told[i][n_nodes_y-1],delta_x);
// // 		// 	}
// // 		// 	if(i==n_nodes_x-1)
// // 		// 	{
// // 		// 		// cen[i][j].FE=Flux_calc(2*Told[i][n_nodes_y-1],2*Tr,delta_x);
// // 		// 		cen[i][j].FE=0.0;
// // 		// 	}
// // 		// 	else
// // 		// 	{
// // 		// 		cen[i][j].FE=Flux_calc(Told[i][n_nodes_y-1],Told[i+1][n_nodes_y-1],delta_x);
// // 		// 	}	
// // 		// 	// TOP WALL TEMPERATURE		
// // 		// 	T[i][n_nodes_y-1] = Told[i][n_nodes_y-1] + (rhs_coefficient/delta_x)*(cen[i][j].FW-cen[i][j].FE)+(rhs_coefficient/delta_y)*(FS-FN)+dt*S[i][n_nodes_y-1];
			
// // 		// }
// // 		// //LEFT and RIGHT
// // 		// for(j=0;j<n_nodes_y;j++)
// // 		// {
// // 		// 	cen[i][j].FW=Flux_calc(2*Tl,2*Told[0][j],delta_x);			
// // 		// 	cen[i][j].FE=Flux_calc(Told[0][j],Told[1][j],delta_x);			
// // 		// 	if(j==0)
// // 		// 	{
// // 		// 		FS=0.0;
// // 		// 	}
// // 		// 	else
// // 		// 	{
// // 		// 		FS=Flux_calc(Told[0][j-1],Told[0][j],delta_y);
// // 		// 	}
// // 		// 	if(j==n_nodes_y-1)
// // 		// 	{
// // 		// 		FN=0.0;
// // 		// 	}
// // 		// 	else
// // 		// 	{
// // 		// 		FN=Flux_calc(Told[0][j],Told[0][j+1],delta_y);
// // 		// 	}			
			
// // 		// 	// LEFT WALL TEMPERATURE
// // 		// 	T[0][j] = Told[0][j] + (rhs_coefficient/delta_x)*(cen[i][j].FW-cen[i][j].FE)+(rhs_coefficient/delta_y)*(FS-FN)+dt*S[0][j];
// // 		// 	// printf("LEFT: FN=%e FS=%e cen[i][j].FE=%e cen[i][j].FW=%e\n Temperature=%e\n",FN,FS,cen[i][j].FE,cen[i][j].FW,T[0][j]);

// // 		// 	// #################################//
// // 		// 	// cen[i][j].FE=Flux_calc(2*Told[n_nodes_x-1][j],2*Tr,delta_x);	
// // 		// 	cen[i][j].FE=0.0;		
// // 		// 	cen[i][j].FW=Flux_calc(Told[n_nodes_x-2][j],Told[n_nodes_x-1][j],delta_x);				;
// // 		// 	if(j==0)
// // 		// 	{
// // 		// 		FS=0.0;
// // 		// 	}
// // 		// 	else
// // 		// 	{
// // 		// 		FS=Flux_calc(Told[n_nodes_x-1][j-1],Told[n_nodes_x-1][j],delta_y);
// // 		// 	}
// // 		// 	if(j==n_nodes_y-1)
// // 		// 	{
// // 		// 		FN=0.0;
// // 		// 	}
// // 		// 	else
// // 		// 	{
// // 		// 		FN=Flux_calc(Told[n_nodes_x-1][j],Told[n_nodes_x-1][j+1],delta_y);
// // 		// 	}			
// // 		// 	// RIGHT WALL TEMPERATURE		
// // 		// 	T[n_nodes_x-1][j] = Told[n_nodes_x-1][j] + (rhs_coefficient/delta_x)*(cen[i][j].FW-cen[i][j].FE)+(rhs_coefficient/delta_y)*(FS-FN)+dt*S[n_nodes_x-1][j];
			
// // 		// }

