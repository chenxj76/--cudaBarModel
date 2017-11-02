// BarModel1.cpp : 定义控制台应用程序的入口点。
//

//#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

//*******可激发系统参数these are parameters about BarModel of exitable system *******
int const nx = 300, ny = 300;//网格点数grid numbers
double const h = 0.02, dx = 0.4, dy = 0.4;//时间、空间步长time and space step 
double const aa = 0.84, bb = 0.07, eps = 0.03, w = 1.5;//可激发系统参数exiatable system parameters

double u[nx + 2][ny + 2], u0[nx + 2][ny + 2], v[nx + 2][ny + 2], v0[nx + 2][ny + 2];
double du1[nx + 2][ny + 2], dv1[nx + 2][ny + 2];
double D;

double sldx, sldy; 
int ncount, i, j, k, m;
FILE *u300,*u301;

double  funcf(double x, double y);
double  funcg(double x, double y);

int main()
{
		sldx = (1.0 / dx)/ dx;
		sldy = (1.0 / dy)/ dy;
		D = 1.0;
		//***init conditions*****	
	for (i = 0; i < nx + 2; i++){
		for (j = 0; j < ny + 2; j++){
			u0[i][j] = 0.0;
			v0[i][j] = 0.0;
		}
	}
	clock_t start, end;
	start = clock();
	for (ncount = 0; ncount < 5001; ncount++){
		printf("%d \n", ncount);
		//****no flux boundary conditions*****
	for (i = 1; i < nx + 1; i++){
		for (j = 1; j < ny + 1; j++){
		u0[0][j] = u0[1][j];
		u0[nx+1][j] = u0[nx][j];
		u0[i][0] = u0[i][1];
		u0[i][ny+1] = u0[i][ny];

		v0[0][j] = v0[1][j];
		v0[nx+1][j] = v0[nx][j];
		v0[i][0] = v0[i][1];
		v0[i][ny+1] = v0[i][ny];
		}
	}
		//*********** Center Differnce for Space *******
		for (i = 1; i < nx + 1; i++){
			for (j = 1; j < ny + 1; j++){
				du1[i][j] = funcf(u0[i][j], v0[i][j]) + D*(sldx*(u0[i + 1][j] + u0[i - 1][j] - 2 * u0[i][j]) + sldy*(u0[i][j+1] + u0[i][j-1] - 2 * u0[i][j]));
				dv1[i][j] = funcg(u0[i][j], v0[i][j]);
			}
		}
		//*****stimulation with a plane  waves****
		if (ncount < 11 ){
			for (i = 1; i < 4; i++){
				for (j = 1; j < ny + 1; j++){
					du1[i][j] = du1[i][j] + 1.5;//*cos(w*ncount*h)+0.5;
				}
			}
		}
		for (i = 1; i < nx + 1; i++){
			for (j = 1; j < ny + 1; j++){
				//***********Forward Euler for Time ******
				u[i][j] = u0[i][j] + h*du1[i][j];
				v[i][j] = v0[i][j] + h*dv1[i][j];
				//***********Update the tow variabls******
				u0[i][j] = u[i][j];
				v0[i][j] = v[i][j];		
			}
		}
		
		//***********trancation 1/2 of the plane wave to generate a spiral wave******
		if (ncount == 1000 ){
			for (i = 1; i < nx + 1; i++){
				for (j = 1; j < (ny / 2)+1; j++){
					u0[i][j] = 0.0;
					v0[i][j] = 0.0;
				}
			}
		}

		if (ncount%50 == 0){
			char filename[100];
			int fileflag=0;
			if (fileflag == 0){
				sprintf(filename, "u%d", ncount);
				u300 = fopen(filename, "w");
				fileflag = 1;
			}
			for (i= 0; i < nx + 2;i++){
				for (j = 0; j < ny + 2;j++){
					fprintf(u300, "%18.7e\t", u0[j][i]);
					if (j == ny+1)  fprintf(u300, "\n");
				}
			}
			if (fileflag == 1){
				fclose(u300);
			}
		}
	}
	end = clock();
	double time_used = (double)(end - start)/CLOCKS_PER_SEC;
	printf("CPUtime=%f", time_used);
//	system("pause");
}

//***********Tow functions are f and g in Partial Defferential equations of the BarModel ******
double funcf(double x, double y)
{
	double funcf_result;
	funcf_result = -x*(x - 1.0)*(x - (y + bb) / aa) / eps;
	return funcf_result;
}
double funcg(double x, double y)
{
	double funcg_result;
	if (x< 1.0/3.0)
		funcg_result = -y;
	else if(x<1.0 || x==1.0)
		funcg_result = 1.0 - 6.75*x*(x - 1.0)*(x - 1.0) - y;
	else
		funcg_result = 1.0 - y;
	return funcg_result;
}
