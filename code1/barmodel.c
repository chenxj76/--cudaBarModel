// BarModel1.c :
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "bar.h"
void File_save(double*h_u0);
int main()
{
	double*h_u,*h_v,*h_u0,*h_v0,*h_du,*h_dv;
	double*d_u,*d_v,*d_u0,*d_v0,*d_du,*d_dv;

	Manage_Memory(1,&h_u,&h_v,&h_u0,&h_v0,&d_u,&d_v,&d_u0,&d_v0,&d_du,&d_dv,&h_du,&h_dv);
	Call_GPU_Init(&d_u0,&d_v0,nx,ny);
	for(int ncount=1; ncount<10; ncount++){
		Call_GPU_Boun(&d_u0,&d_v0);		
		Call_GPU_Space(&d_u0,&d_v0,&d_du,&d_dv);
		if(ncount<10)Call_GPU_PlaneWave(&d_du,ncount);		
		Call_GPU_ForEuler(&d_u,&d_u0,&d_du,&d_v,&d_v0,&d_dv);//调用这句的时候du[idx]数据为什么已经清零了？
		Call_GPU_Update(&d_u,&d_u0,&d_v,&d_v0);
		if(ncount==1000)Call_GPU_Trancation(&d_u0,&d_v0);
		if(ncount==5)Manage_Comms(2,&h_u0,&d_u0);
	}	
	File_save(h_u0);
	Manage_Memory(2,&h_u,&h_v,&h_u0,&h_v0,&d_u,&d_v,&d_u0,&d_v0,&d_du,&d_dv,&h_du,&h_dv);
	return 0;
}

//--------------------------------------------------------------------------
void File_save(double*h_u0){
//输出文件保存代码块
	FILE *pFile1;
	pFile1 = fopen("Upattern1", "w");
	int idx;
	for (idx= 0; idx <(nx+1)*(ny+1); idx++){
	//	printf("u[%d]=%18.7e",idx,h_u0[idx]);
		fprintf(pFile1, "%18.7e", h_u0[idx]);		
	if ((idx+1)%nx==0)  fprintf(pFile1, "\n");
	}				
	fclose(pFile1);
	//------------------------------------------
	FILE *pFile;
	pFile = fopen("Upattern", "w");
	int ix,iy;
	double u0[301][301];	
	for (int i= 0; i <nx; i++){
	for (int j= 0; j <ny; j++){
		ix=i;
		iy=j+1;
		idx = ix*nx+iy;	
		u0[i][j]=h_u0[idx];  					//数组的角标容易出错，globalidx的xy与数组的index的xy刚好相反
		//printf("u[i][j]=%18.7e",u0[i][j]);
		fprintf(pFile, "%18.7e", u0[i][j]);		
	}//printf("\n");
	if (i!=(nx-1))  fprintf(pFile, "\n");
	}				
	fclose(pFile);
	//---------------------------------------------
}
