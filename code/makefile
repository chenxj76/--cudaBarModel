CUDA_INSTALL_PATH:=/usr/local/cuda
CUDA_LIB:=/usr/local/cuda/lib64 -lcuda -lcudart

all:CPU GPU1 GPU2
	g++ barmodel.o bar.o memory.o -o3 -o bar.run -L $(CUDA_LIB)
CPU:
	g++ barmodel.c -c
GPU1:
	nvcc -arch=sm_35 bar.cu -c
GPU2:
	nvcc -arch=sm_35 memory.cu -c
clean:
	rm *.o