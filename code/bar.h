#define BLOCK_SIZE 32

//*******these are parameters about BarModel of exitable system *******
#define nx 300
#define ny 300
#define h 0.02//time step 
#define dx 0.4//space step 
#define dy 0.4
#define aa 0.84
#define bb 0.07
#define eps 0.03
#define w 1.5
#define D 1.0

void Manage_Memory(int phase);
void Manage_Comms(int phase);
void Call_GPU_Init();
void Call_GPU_Boun();
void Call_GPU_Space();
void Call_GPU_PlaneWave(int ncount);
void Call_GPU_ForEuler();
void Call_GPU_Update();
void Call_GPU_Trancation();