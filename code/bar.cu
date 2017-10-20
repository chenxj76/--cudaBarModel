#include "bar.h"
#include <cuda_runtime.h>
#include <stdio.h>
__global__ void GPUInitOnDevice(double*u,double*v,int x,int y);
__global__ void GPUCalcBar(double*u,double*v,double*u0,double*v0,double*du,double*dv);

void Manage_Memory(int phase,double**u,double**v,double**u0,double**v0,double**dd_u,double**dd_v,double**dd_u0,double**dd_v0,double**dd_du,double**dd_dv){
	cudaError_t Error;
	size_t size=nx*ny*sizeof(double);
	if (phase==1){
		*u=(double*)malloc(size);
		*v=(double*)malloc(size);
		*u0=(double*)malloc(size);
		*v0=(double*)malloc(size);
		Error=cudaMalloc((void**)dd_u,size);
		Error=cudaMalloc((void**)dd_u0,size);
		Error=cudaMalloc((void**)dd_v,size);
		Error=cudaMalloc((void**)dd_v0,size);
		Error=cudaMalloc((void**)dd_du,size);
		Error=cudaMalloc((void**)dd_dv,size);
		printf("MemoryMalloc:%s\n",cudaGetErrorString(Error));
	}
	if (phase==2){
		free(*u);
		free(*v);
		free(*u0);
		free(*v0);
		Error=cudaFree(*dd_u);
		Error=cudaFree(*dd_u0);
		Error=cudaFree(*dd_v);
		Error=cudaFree(*dd_v0);
		Error=cudaFree(*dd_du);
		Error=cudaFree(*dd_dv);
	}

}

void Manage_Comms(int phase,double*hh_u,double*dd_u){

	cudaError_t Error;
	size_t size=nx*ny*sizeof(double);
	if (phase==2){				
	cudaMemcpy(hh_u,dd_u,size,cudaMemcpyDeviceToHost);
		Error=cudaMemcpy(hh_u,dd_u,size,cudaMemcpyDeviceToHost);
		printf("device to host:%s\n",cudaGetErrorString(Error));
		}
	if (phase==1){
		Error=cudaMemcpy(dd_u,hh_u,size,cudaMemcpyHostToDevice);
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
void Call_GPU_Init(double*u0,double*v0,int x,int y){	
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUInitOnDevice<<<grid,block>>>(u0,v0,x,y);
	Error=cudaDeviceSynchronize();
	printf("InitDeviceSynchronize:%s\n",cudaGetErrorString(Error));
}

__global__ void GPUCalcBar(double*u,double*v,double*u0,double*v0,double*du,double*dv){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; //matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx
unsigned int idx = iy*nx+ ix;//globalIdx
		double*hh_u0;
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
		if (ncount == 5000){
		 *hh_u0=*u0;
		}
	
	}	
}
void Call_GPU_Calc_Bar(double*u,double*v,double*u0,double*v0,double*du,double*dv){
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUCalcBar<<<grid,block>>>(u,v,u0,v0,du,dv);
	Error=cudaDeviceSynchronize();
	printf("CalcDeviceSynchronize:%s\n",cudaGetErrorString(Error));
}




