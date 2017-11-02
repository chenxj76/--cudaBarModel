// BarModel1.c :
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "bar.h"
#include "systime.h"
extern double*h_u,*h_v,*h_u0,*h_v0,*h_du,*h_dv;
extern double*d_u,*d_v,*d_u0,*d_v0,*d_du,*d_dv;
void File_save(int n);
double Start,End;
int main()
{
	Manage_Memory(1);	
	Call_GPU_Init();//直接在GPU里初始化，会有问题吗？
	Start=cpuSecond();
	for(int ncount=0;ncount<5001;ncount++){
		Call_GPU_Boun();
		Call_GPU_Space();
		if(ncount<11)Call_GPU_PlaneWave(ncount);
		Call_GPU_ForEuler();
		Call_GPU_Update();		
		if(ncount==1000)Call_GPU_Trancation();
		if(ncount%50==0){
			Manage_Comms(2);
			int fileflag = 0;
			File_save(ncount);
		}
	}
	//Manage_Comms(2);
	End=cpuSecond()-Start;//这样子计算GPU的时间是不是不准确，因为包括了数据从GPU送到CPU的时间。
	printf("GPUtime=%f\n",End);
	//File_save();
	Manage_Memory(2);
	return 0;
}

void File_save(int n){

	FILE *ap;
	char filename[100];
	int ix,iy,idx;
	int fileflag;
	if (fileflag == 0){
						sprintf(filename, "ap%d", n);
						ap = fopen(filename, "w");
						fileflag = 1;
					}
	for (idx= 0; idx <(nx+2)*(ny+2); idx++){
		fprintf(ap, "%18.7e\t", h_u0[idx]);		
		if ((idx+1)%(nx+2)==0)  fprintf(ap, "\n");
	}	
	if (fileflag == 1){
			fclose(ap);
		}	
	/*	
	for (int i= 0; i <nx+2; i++){
		for (int j= 0; j <ny+2; j++){ 
			int id=i*(nx+2)+j;//扩展板的globalIdx,
			idx = (i+1)*nx+ (j+1)+2*i;//globalIdx,不要计算4条边界	
			u0[i][j]=h_u0[idx];  	//数组的角标容易出错，globalidx的xy与数组的index的xy刚好相反
			//printf("u[i][j]=%18.7e",u0[i][j]);
			//if(i!=0&&j!=0&&i!=nx+1&&j!=ny+1){
				fprintf(pFile, "%18.7e", u0[i][j]);
			if (j!=ny+1) fprintf(pFile, "\t");
			else 	   fprintf(pFile, "\n");
			
		}
			//printf("\n");
							
	}	*/	
	
}
