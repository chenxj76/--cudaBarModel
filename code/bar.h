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

void Manage_Memory(int phase,double**u,double**v,double**u0,double**v0,double**dd_u,double**dd_v,double**dd_u0,double**dd_v0,double**dd_du,double**dd_dv);
void Manage_Comms(int phase,double*hh_u,double*dd_u);
void Call_GPU_Init(double*u0,double*v0,int x,int y);
double Func_FG(int phase,double*x,double*y);
void Call_GPU_Calc_Bar(double*u,double*v,double*u0,double*v0,double*du,double*dv);
//void Save_Results(double*hh_u0);