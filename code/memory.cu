#include "bar.h"
#include <cuda_runtime.h>
#include <stdio.h>

double*h_u,*h_v,*h_u0,*h_v0,*h_du,*h_dv;
double*d_u,*d_v,*d_u0,*d_v0,*d_du,*d_dv;


void Manage_Memory(int phase){
	cudaError_t Error;
	size_t size=(nx+2)*(ny+2)*sizeof(double);
	if (phase==1){
		h_u=(double*)malloc(size);
		h_v=(double*)malloc(size);
		h_u0=(double*)malloc(size);
		h_v0=(double*)malloc(size);
		Error=cudaMalloc((void**)&d_u,size);
		Error=cudaMalloc((void**)&d_u0,size);
		Error=cudaMalloc((void**)&d_v,size);
		Error=cudaMalloc((void**)&d_v0,size);
		Error=cudaMalloc((void**)&d_du,size);
		Error=cudaMalloc((void**)&d_dv,size);
		printf("MemoryMalloc:%s\n",cudaGetErrorString(Error));
	}
	if (phase==2){
		free(h_u);
		free(h_v);
		free(h_u0);
		free(h_v0);
		Error=cudaFree(d_u);
		Error=cudaFree(d_u0);
		Error=cudaFree(d_v);
		Error=cudaFree(d_v0);
		Error=cudaFree(d_du);
		Error=cudaFree(d_dv);
	}

}

void Manage_Comms(int phase){

	cudaError_t Error;
	size_t size=(nx+2)*(ny+2)*sizeof(double);
	if (phase==2){				
		Error=cudaMemcpy(h_u0,d_u0,size,cudaMemcpyDeviceToHost);
		if (Error != cudaSuccess)printf("device to host:%s\n",cudaGetErrorString(Error));
		}
	if (phase==1){
		Error=cudaMemcpy(d_u0,h_u0,size,cudaMemcpyHostToDevice);
		if (Error != cudaSuccess)printf("host to device:%s\n",cudaGetErrorString(Error));}
}
