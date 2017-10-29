#define BLOCK_SIZE 30
#define TILE_SIZE 10
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

void Manage_Memory(int phase,double**h_u,double**h_v,double**h_u0,double**h_v0,double**d_u,double**d_v,double**d_u0,double**d_v0,double**d_du,double**d_dv,double**h_du,double**h_dv);
void Manage_Comms(int phase,double**h_u0,double**d_u0);
void Call_GPU_Init(double**d_u0,double**d_v0,int x,int y);
void Call_GPU_Boun(double**d_u0,double**d_v0);
void Call_GPU_Space(double**d_u0,double**d_v0,double**d_du,double**d_dv);
void Call_GPU_PlaneWave(double**d_du,int d_ncount);
void Call_GPU_ForEuler(double**d_u,double**d_u0,double**d_du,double**d_v,double**d_v0,double**d_dv);
void Call_GPU_Update(double**d_u,double**d_u0,double**d_v,double**d_v0);
void Call_GPU_Trancation(double**d_u0,double**d_v0);