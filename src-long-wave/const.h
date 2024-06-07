/*
Programme to solve the 1D Long-Wave OWT equation

i PSI_t = - (1/2) PSI_xx    -  |PSI|^2 PSI  -  (1/g)* PSI |PSI|^2_xx      + dissipation + forcing

Author: Jason Laurie
Date: 31/01/2022

*/

const int N = 4096;                      //number of spatial grid points
const int M = 3*N/2;                      //defined the size of the anti-aliasing array M=3N/2
const double g = 512.0;
const int num_threads = 1;


//===============set timestepping parameters============================

const bool FLAG_TIMESTEP_ETDRK = true;  //if true routine is ETDRK4 else RK2
const int FLAG_ETDRK_ORDER = 4;
const double dt = 1.e-7;                //time step
const int nsteps = 999999999;           //number of total time steps
const int outstep = 5e4;                //outputs data at these time steps
//======================parameters for domain ===============================
const double pi = 3.14159265358979323846;   //pi
const double Lx = 2.0*pi;                   //length of the box
const double dx = Lx/ double(N);            //grid size


//===========parameters from equation =============================
const bool FLAG_HYPER_DISS = true;
const double nu = 1.e-38;//1.e-32;					//hyperviscosity coefficient
const double nupower = 6.0;					//power of the laplacian 	

const bool FLAG_HYPO_DISS = false;
const double alpha = 0.e2;					//friction coefficient
const double alphapower =-2.0;				//negative power of the laplacian


//===================parameters for forcing routine ===============================
const bool FLAG_FORCE_ON = false; //Flag to turn on additive forcing
const bool FLAG_FORCE_AMP_RESCALE = false; // Flag for having forcing amplitude rescaled to give unit energy density
const bool FLAG_FORCE_EXP = false;  // turn to true to have exponential forcing spectrum **default is annulus***
const double kf = 256;							//forcing wavenumber
const double dk = 4.0;							//width of forcing annulus (in terms of number of modes)
const double Amp = 1.e-0;//0.001;						//amplitude of forcing
const double force_power = 4.0; // for exp force only f_k = (k/k_f)^force_power * exp( -(k / k_f)^force_power )
const double mean_waveaction = 1.0; //paramter for rescaling forcing wih large-scale dissipation

