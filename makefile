CUDA_INSTALL_PATH:=/usr/local/cuda
CUDA_LIB:=/usr/local/cuda/lib64 -lcuda -lcudart

all:CPU GPU
	g++ barmodel.o bar.o -o3 -o bar.run -L $(CUDA_LIB)
CPU:
	g++ barmodel.c -c
GPU:
	nvcc -arch=sm_35 bar.cu -c
clean:
<<<<<<< HEAD
	rm *.o
=======
	rm *.o
>>>>>>> d23285bd597c8ba23bf4d16cb949eb997b15c9b7
