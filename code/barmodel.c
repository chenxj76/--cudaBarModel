// BarModel1.c :
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "bar.h"

int main()
{
	double*h_u,*h_v,*h_u0,*h_v0,*h_du,*h_dv;
	double*d_u,*d_v,*d_u0,*d_v0,*d_du,*d_dv;
	Manage_Memory(1,&h_u,&h_v,&h_u0,&h_v0,&d_u,&d_v,&d_u0,&d_v0,&d_du,&d_dv);
	Call_GPU_Init(&d_u0,&d_v0,nx,ny);
	Call_GPU_Calc_Bar(&d_u,&d_v,&d_u0,&d_v0,&d_du,&d_dv);
	Manage_Comms(2,&h_u0,&d_u0);
	//输出文件保存代码块
	FILE *pFile;
	pFile = fopen("Upattern", "w");
	for (int idx= 0; idx < nx*ny; idx++){
      if (idx%nx==0 && idx>0)  fprintf(pFile, "\n");
	  fprintf(pFile, "%18.7e", h_u0[idx]);	
	}
	fclose(pFile);
	Manage_Memory(2,&h_u,&h_v,&h_u0,&h_v0,&d_u,&d_v,&d_u0,&d_v0,&d_du,&d_dv);
	return 0;
}


