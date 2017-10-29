#include "bar.h"
#include <cuda_runtime.h>
#include <stdio.h>
double*h_u,*h_v,*h_u0,*h_v0,*h_du,*h_dv;
double*d_u,*d_v,*d_u0,*d_v0,*d_du,*d_dv;

__global__ void GPUInitOnDevice(double*u,double*v,int x,int y);
__global__ void GPUCalcBar(double*u,double*v,double*u0,double*v0,double*du,double*dv);
__global__ void GPUBoundary(double*u,double*v);
__global__ void GPUCalcSpace(double*u0,double*v0,double*du,double*dv);
__global__ void GPUPlaneWave(double*du,int ncount);
__global__ void GPUForEuler(double*u,double*u0,double*du,double*v,double*v0,double*dv);
__global__ void GPUUpdate(double*u,double*u0,double*v,double*v0);
__global__ void GPUTrancation(double*u,double*v);

void Manage_Memory(int phase,double**h_u,double**h_v,double**h_u0,double**h_v0,double**d_u,double**d_v,double**d_u0,double**d_v0,double**d_du,double**d_dv,double**h_du,double**h_dv){
	cudaError_t Error;
	size_t size=(nx+1)*(ny+1)*sizeof(double);//index需要+1，因为使用了no flux 边界，需要往外扩大一条边。
	if (phase==1){
		*h_u=(double*)malloc(size);
		*h_v=(double*)malloc(size);
		*h_u0=(double*)malloc(size);
		*h_v0=(double*)malloc(size);
		*h_du=(double*)malloc(size);
		*h_dv=(double*)malloc(size);
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
		free(*h_du);
		free(*h_dv);
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
	size_t size=(nx+1)*(ny+1)*sizeof(double);
	if (phase==2){				
		Error=cudaMemcpy(*h_u0,*d_u0,size,cudaMemcpyDeviceToHost);
		printf("device to host:%s\n",cudaGetErrorString(Error));
		}
	if (phase==1){
		Error=cudaMemcpy(*d_u0,*h_u0,size,cudaMemcpyHostToDevice);
		printf("host to device:%s\n",cudaGetErrorString(Error));}
}
//----------------Init-------------------
__global__ void GPUInitOnDevice(double*u,double*v,int x,int y){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;//matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx
unsigned int idx = iy*x+ ix;//globalIdx
	//***init conditions*****			
		if(ix<(x+1)&&iy<(y+1)){
			u[idx] = 1.0; 
			v[idx] = 1.0;			
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
//----------------boundary-------------------
__global__ void GPUBoundary(double*u,double*v){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; //matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx
		//****no flux boundary conditions*****
				u[iy*nx] = u[1+iy*nx];
				u[nx+iy*nx] = u[(nx-1)+iy*nx];
				v[iy*nx] = v[1+iy*nx];
				v[nx+iy*nx] = v[(nx-1)+iy*nx];
				u[ix] = u[ix+nx];
				u[ix+ny*nx] = u[ix+(ny-1)*nx];
				v[ix] = v[ix+nx];
				v[ix+ny*nx] = v[ix+(ny-1)*nx];		
}		
void Call_GPU_Boun(double**d_u0,double**d_v0){	
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUBoundary<<<grid,block>>>(*d_u0,*d_v0);
	Error=cudaDeviceSynchronize();
	printf("GPUBoundarySynchronize:%s\n",cudaGetErrorString(Error));
}			
			
//----------------Center Differnce for Space-------------------			
__global__ void GPUCalcSpace(double*u0,double*v0,double*du,double*dv){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; //matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx
unsigned int idx = iy*nx+ ix;//globalIdx
		double sldx, sldy; 
		sldx = (1.0 / dx)/ dx;
		sldy = (1.0 / dy)/ dy;
		
		//if(ix>0&&iy>0){//此处需要控制边界
			du[idx] = -u0[idx]*(u0[idx] - 1.0)*(u0[idx] - (v0[idx] + bb) / aa) / eps + D*(sldx*(u0[(ix+1)+iy*nx] + u0[(ix-1)+iy*nx] - 2 * u0[ix+iy*nx]) 
									+ sldy*(u0[ix+(iy+1)*nx] + u0[ix+(iy-1)*nx] - 2 * u0[ix+iy*nx]));
			
			{
			if(u0[idx]< 1.0/3.0) dv[idx] = -v0[idx];
			else if(u0[idx]<1.0 || u0[idx]==1.0)
			dv[idx] = 1.0 - 6.75*u0[idx]*(u0[idx] - 1.0)*(u0[idx] - 1.0) - v0[idx];
			else dv[idx] = 1.0 - v0[idx];
						
			}	
		//}
}		
void Call_GPU_Space(double**d_u0,double**d_v0,double**d_du,double**d_dv){	
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUCalcSpace<<<grid,block>>>(*d_u0,*d_v0,*d_du,*d_dv);
	Error=cudaDeviceSynchronize();
	printf("GPUCalcSpaceSynchronize:%s\n",cudaGetErrorString(Error));
}			

//----------------stimulation with a plane waves-------------------		
__global__ void GPUPlaneWave(double*du,int ncount){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; //matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx
unsigned int idx = iy*nx+ ix;//globalIdx	
		if (ncount < 10 && ix<3 ){
				du[idx] = du[idx] + 1.5*cos(w*ncount*h) + 0.5;
				
		}
}
void Call_GPU_PlaneWave(double**d_du,int ncount){	
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUPlaneWave<<<grid,block>>>(*d_du,ncount);
	Error=cudaDeviceSynchronize();
	printf("GPUPlaneWaveSynchronize:%s\n",cudaGetErrorString(Error));
}
//----------------Forward Euler for Time-------------------
__global__ void GPUForEuler(double*u,double*u0,double*du,double*v,double*v0,double*dv){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; //matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx
unsigned int idx = iy*nx+ ix;//globalIdx
				printf("%18.7e",du[idx]);		
				u[idx] = u0[idx] + h*du[idx];
				v[idx] = v0[idx] + h*dv[idx];
}	
void Call_GPU_ForEuler(double**d_u,double**d_u0,double**d_du,double**d_v,double**d_v0,double**d_dv){	
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUForEuler<<<grid,block>>>(*d_u,*d_u0,*d_du,*d_v,*d_v0,*d_dv);
	Error=cudaDeviceSynchronize();
	printf("GPUForEulerSynchronize:%s\n",cudaGetErrorString(Error));
}			
//----------------Update the tow variabls-------------------
__global__ void GPUUpdate(double*u,double*u0,double*v,double*v0){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; //matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx
unsigned int idx = iy*nx+ ix;//globalIdx	
				u0[idx] = u[idx];
				v0[idx] = v[idx];
}
void Call_GPU_Update(double**d_u,double**d_u0,double**d_v,double**d_v0){	
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUUpdate<<<grid,block>>>(*d_u,*d_u,*d_v,*d_v0);
	Error=cudaDeviceSynchronize();
	printf("GPUUpdateSynchronize:%s\n",cudaGetErrorString(Error));
}				
//----------------trancation 1/2 of the plane wave to generate a spiral wave-------------------
__global__ void GPUTrancation(double*u,double*v){
unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; //matrixIdx
unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;//matrixIdx
unsigned int idx = iy*nx+ ix;//globalIdx			
		if (idx<nx+1 && iy<(ny/2)+1){
				u[idx] = 0;
				v[idx] = 0;
		}
}
void Call_GPU_Trancation(double**d_u0,double**d_v0){	
	dim3 block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 grid((nx+block.x-1)/block.x,(ny+block.y-1)/block.y);
	cudaError_t Error;
	GPUTrancation<<<grid,block>>>(*d_u0,*d_v0);
	Error=cudaDeviceSynchronize();
	printf("GPUTrancationSynchronize:%s\n",cudaGetErrorString(Error));
}



