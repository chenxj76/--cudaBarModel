#include "bar.h"
#include <cuda_runtime.h>
#include <stdio.h>
double*h_u,*h_v,*h_u0,*h_v0,*h_du,*h_dv;
double*d_u,*d_v,*d_u0,*d_v0,*d_du,*d_dv;
__global__ void GPUInitOnDevice(double*u,double*v,int x,int y);
__global__ void GPUCalcBar(double*u,double*v,double*u0,double*v0,double*du,double*dv);

void Manage_Memory(int phase,double**h_u,double**h_v,double**h_u0,double**h_v0,double**d_u,double**d_v,double**d_u0,double**d_v0,double**d_du,double**d_dv){
	cudaError_t Error;
	size_t size=nx*ny*sizeof(double);
	if (phase==1){
		*h_u=(double*)malloc(size);
		*h_v=(double*)malloc(size);
		*h_u0=(double*)malloc(size);
		*h_v0=(double*)malloc(size);
		Error=cudaMalloc((void**)d_u,size);
		Error=cudaMalloc((void**)d_u0,size);
		Error=cudaMalloc((void**)d_v,size);
		Error=cudaMalloc((void**)d_v0,size);
		Error=cudaMalloc((void**)d_du,size);
		Error=cudaMalloc((void**)d_dv,size);
		printf("MemoryMalloc:%s\n",cudaGetErrorString(Error));
	}
	if (phase==2){
		free(*h_u);
		free(*h_v);
		free(*h_u0);
		free(*h_v0);
		Error=cudaFree(*d_u);
		Error=cudaFree(*d_u0);
		Error=cudaFree(*d_v);
		Error=cudaFree(*d_v0);
		Error=cudaFree(*d_du);
		Error=cudaFree(*d_dv);
	}

}

void Manage_Comms(int phase,double**h_u0,double**d_u0){

	cudaError_t Error;
	size_t size=nx*ny*sizeof(double);
	if (phase==2){				
		Error=cudaMemcpy(*h_u0,*d_u0,size,cudaMemcpyDeviceToHost);
		printf("device to host:%s\n",cudaGetErrorString(Error));
		}
	if (phase==1){
		Error=cudaMemcpy(*d_u0,*h_u0,size,cudaMemcpyHostToDevice);
		printf("host to device:%s\n",cudaGetErrorString(Error));}
}

__global__ void GPUInitOnDevice(double*u,double*v,int x,int y){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;//matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx
unsigned int idx = iy*x+ ix;//globalIdx
	//***init conditions*****			
		if(ix<x&&iy<y){
			u[idx] = 0.0; 
			v[idx] = 0.0;			
	}
}
void Call_GPU_Init(double**d_u0,double**d_v0,int x,int y){	
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUInitOnDevice<<<grid,block>>>(*d_u0,*d_v0,x,y);
	Error=cudaDeviceSynchronize();
	printf("InitDeviceSynchronize:%s\n",cudaGetErrorString(Error));
}

__global__ void GPUCalcBar(double*u,double*v,double*u0,double*v0,double*du,double*dv){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; //matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx
unsigned int idx = iy*nx+ ix;//globalIdx

		double sldx, sldy; 
		int ncount;
		sldx = (1.0 / dx)/ dx;
		sldy = (1.0 / dy)/ dy;

	for (ncount = 1; ncount < 10001; ncount++){

		//****no flux boundary conditions*****
				u0[iy*nx] = u0[1+iy*nx];
				u0[nx+iy*nx] = u0[(nx-1)+iy*nx];
				v0[iy*nx] = v0[1+iy*nx];
				v0[nx+iy*nx] = v0[(nx-1)+iy*nx];
				u0[ix] = u0[ix+nx];
				u0[ix+ny*nx] = u0[ix+(ny-1)*nx];
				v0[ix] = v0[ix+nx];
				v0[ix+ny*nx] = v0[ix+(ny-1)*nx];		
				
		//*********** Center Differnce for Space *******
		
		
			du[idx] = -u0[idx]*(u0[idx] - 1.0)*(u0[idx] - (v0[idx] + bb) / aa) / eps + D*(sldx*(u0[(ix+1)+iy*nx] + u0[(ix-1)+iy*nx] - 2 * u0[ix+iy*nx]) 
									+ sldy*(u0[ix+(iy+1)*nx] + u0[ix+(iy-1)*nx] - 2 * u0[ix+iy*nx]));
			
			{
			if (u0[idx]< 1.0/3.0)
			dv[idx] = -v0[idx];
			else if(u0[idx]<1.0 || u0[idx]==1.0)
			dv[idx] = 1.0 - 6.75*u0[idx]*(u0[idx] - 1.0)*(u0[idx] - 1.0) - v0[idx];
			else
			dv[idx] = 1.0 - v0[idx];
			}

			
		//*****stimulation with a plane waves****
		if (ncount < 10 && ix<3 ){
				du[idx] = du[idx] + 1.5*cos(w*ncount*h) + 0.5;
		}

		//***********Forward Euler for Time ******
				u[idx] = u0[idx] + h*du[idx];
				v[idx] = v0[idx] + h*dv[idx];
		//***********Update the tow variabls******
				u0[idx] = u[idx];
				v0[idx] = v[idx];		
		
		//***********trancation 1/2 of the plane wave to generate a spiral wave******
		if (ncount == 1000 && idx<nx+1 && iy<(ny/2)+1){
				u0[idx] = 0;
				v0[idx] = 0;
		}
	//	if (ncount == 5000){
	//	 h_u0[idx]=u0[idx];   //问题出在这CalcDeviceSynchronize:an illegal memory access was encountered！
	//	}
	
	}	
}
void Call_GPU_Calc_Bar(double**d_u,double**d_v,double**d_u0,double**d_v0,double**d_du,double**d_dv){
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUCalcBar<<<grid,block>>>(*d_u,*d_v,*d_u0,*d_v0,*d_du,*d_dv);
	Error=cudaDeviceSynchronize();
	printf("CalcDeviceSynchronize:%s\n",cudaGetErrorString(Error));
}




