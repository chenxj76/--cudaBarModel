#include "bar.h"
#include <cuda_runtime.h>
#include <stdio.h>
extern double*h_u,*h_v,*h_u0,*h_v0,*h_du,*h_dv;
extern double*d_u,*d_v,*d_u0,*d_v0,*d_du,*d_dv;


__global__ void GPUInitOnDevice(int x,int y,double*d_u,double*d_v);
__global__ void GPUBoundary(double*d_u,double*d_v);
__global__ void GPUCalcSpace(double*d_u0,double*d_v0,double*d_du,double*d_dv);
__global__ void GPUPlaneWave(double*d_du,int ncount);
__global__ void GPUForEuler(double*d_u,double*d_u0,double*d_du,double*d_v,double*d_v0,double*d_dv);
__global__ void GPUUpdate(double*d_u,double*d_u0,double*d_v,double*d_v0);
__global__ void GPUTrancation(double*d_u,double*d_v);


//----------------Init-------------------
__global__ void GPUInitOnDevice(double*d_u,double*d_v){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;//matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx

	//***init conditions*****			
		if(ix<nx+2&&iy<nx+2){//条件控制有多少个点赋初值
		unsigned int idx = iy*(nx+2)+ ix;//globalIdx	
			d_u[idx] = 0.0; 
			d_v[idx] = 0.0;			
	}
}
void Call_GPU_Init(){	
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+2+block.x-1)/block.x,(ny+2+block.y-1)/block.y);
	cudaError_t Error;
	GPUInitOnDevice<<<grid,block>>>(d_u,d_v);
	Error=cudaDeviceSynchronize();
	if (Error != cudaSuccess)printf("InitDeviceSynchronize:%s\n",cudaGetErrorString(Error));
}
//----------------boundary-------------------
__global__ void GPUBoundary(double*d_u,double*d_v){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; //matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx
//unsigned int id = iy*(nx+2)+ ix;//globalIdx，包括4条边
		
		//****no flux boundary conditions*****
		if(ix<nx&&iy<nx){
				d_u[(iy+1)*(nx+2)+0] = d_u[(iy+1)*(nx+2)+1];
				d_u[(iy+1)*(nx+2)+(nx+1)] = d_u[(iy+1)*(nx+2)+nx];
				d_u[(ix+1)+0*(nx+2)] = d_u[(ix+1)+1*(nx+2)];
				d_u[(ix+1)+(ny+1)*(nx+2)] = d_u[(ix+1)+ny*(nx+2)];
				
				d_v[(iy+1)*(nx+2)+0] = d_v[(iy+1)*(nx+2)+1];
				d_v[(iy+1)*(nx+2)+(nx+1)] = d_v[(iy+1)*(nx+2)+nx];			
				d_v[(ix+1)+0*(nx+2)] = d_v[(ix+1)+1*(nx+2)];
				d_v[(ix+1)+(ny+1)*(nx+2)] = d_v[(ix+1)+ny*(nx+2)];
		}									
}		
void Call_GPU_Boun(){	
	dim3 block(BLOCK_SIZE,1);
	dim3 grid((nx+block.x-1)/block.x,1);
	cudaError_t Error;
	GPUBoundary<<<grid,block>>>(d_u,d_v);
	Error=cudaDeviceSynchronize();
	if (Error != cudaSuccess)printf("GPUBoundarySynchronize:%s\n",cudaGetErrorString(Error));
}			
			
//----------------Center Differnce for Space-------------------			
__global__ void GPUCalcSpace(double*d_u0,double*d_v0,double*d_du,double*d_dv){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; //matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx
//unsigned int id=iy*(nx+2)+ix;//扩充了的global index
		double sldx, sldy; 
		sldx = (1.0 / dx)/ dx;
		sldy = (1.0 / dy)/ dy;

		if(ix<nx&&iy<ny){
			
			int idx = (iy+1)*(nx+2)+ (ix +1);//globalIdx,不要计算4条边界
			d_du[idx] = -d_u0[idx]*(d_u0[idx] - 1.0)*(d_u0[idx] - (d_v0[idx] + bb) / aa) / eps + D*(sldx*(d_u0[idx+1] + d_u0[idx-1] - 2 * d_u0[idx]) 
									+ sldy*(d_u0[idx+nx+2] + d_u0[idx-nx-2] - 2 * d_u0[idx]));
			
			
			if(d_u0[idx]< 1.0/3.0) d_dv[idx] = -d_v0[idx];
			else if(d_u0[idx]<1.0 || d_u0[idx]==1.0)
			d_dv[idx] = 1.0 - 6.75*d_u0[idx]*(d_u0[idx] - 1.0)*(d_u0[idx] - 1.0) - d_v0[idx];
			else d_dv[idx] = 1.0 - d_v0[idx];
						
				
		}
}		
void Call_GPU_Space(){	
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUCalcSpace<<<grid,block>>>(d_u0,d_v0,d_du,d_dv);
	Error=cudaDeviceSynchronize();
	if (Error != cudaSuccess)printf("GPUCalcSpaceSynchronize:%s\n",cudaGetErrorString(Error));
}			

//----------------stimulation with a plane waves-------------------		
__global__ void GPUPlaneWave(double*d_du,int ncount){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; //matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx

	
		if (iy < nx && ix<3 ){
				
				int idx = (iy+1)*(nx+2)+ (ix +1);//globalIdx,不要计算4条边界
				d_du[idx] = d_du[idx] + 1.5*cos(w*ncount*h) + 0.5;
				
		}
}
void Call_GPU_PlaneWave(int ncount){	
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUPlaneWave<<<grid,block>>>(d_du,ncount);
	Error=cudaDeviceSynchronize();
	if (Error != cudaSuccess)printf("GPUPlaneWaveSynchronize:%s\n",cudaGetErrorString(Error));
}
//----------------Forward Euler for Time-------------------
__global__ void GPUForEuler(double*d_u,double*d_u0,double*d_du,double*d_v,double*d_v0,double*d_dv){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; //matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx

		if(ix<nx&&iy<ny){
		
				int idx = (iy+1)*(nx+2)+ (ix +1);//globalIdx,不要计算4条边界	
				d_u[idx] = d_u0[idx] + h*d_du[idx];
				d_v[idx] = d_v0[idx] + h*d_dv[idx];
				}
}	
void Call_GPU_ForEuler(){	
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUForEuler<<<grid,block>>>(d_u,d_u0,d_du,d_v,d_v0,d_dv);
	Error=cudaDeviceSynchronize();
	if (Error != cudaSuccess)printf("GPUForEulerSynchronize:%s\n",cudaGetErrorString(Error));
}			
//----------------Update the tow variabls-------------------
__global__ void GPUUpdate(double*d_u,double*d_u0,double*d_v,double*d_v0){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; //matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx

		if(ix<nx&&iy<ny){
				
				int idx = (iy+1)*(nx+2)+ (ix +1);//globalIdx,不要计算4条边界
				d_u0[idx] = d_u[idx];
				d_v0[idx] = d_v[idx];
				}
}
void Call_GPU_Update(){	
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUUpdate<<<grid,block>>>(d_u,d_u0,d_v,d_v0);
	Error=cudaDeviceSynchronize();
	if (Error != cudaSuccess)printf("GPUUpdateSynchronize:%s\n",cudaGetErrorString(Error));
}				
//----------------trancation 1/2 of the plane wave to generate a spiral wave-------------------
__global__ void GPUTrancation(double*d_u0,double*d_v0){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; //matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx
unsigned int idx;//globalIdx
		if(ix<nx&&iy<(ny/2)+1){
				idx = (iy+1)*(nx+2)+ (ix +1);//globalIdx,不要计算4条边界			
				d_u0[idx] = 0.0;
				d_v0[idx] = 0.0;
		}
}
void Call_GPU_Trancation(){	
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUTrancation<<<grid,block>>>(d_u0,d_v0);
	Error=cudaDeviceSynchronize();
	if (Error != cudaSuccess)printf("GPUTrancationSynchronize:%s\n",cudaGetErrorString(Error));
}








