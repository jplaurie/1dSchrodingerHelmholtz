/*
Programme to solve the 1D Schrödinger-Helmholtz / nonlinear Schrödinger equation

"SHE"

i PSI_t = c * PSI_xx    +  g * PSI (1 - beta * d^2_xx)^-1 |PSI|^2  + mu * PSI    + dissipation + forcing

or

"NLS"

i PSI_t = c * PSI_xx    +  g * PSI |PSI|^2  +   mu * PSI    + dissipation + forcing

Author: Jason Laurie
Date: 27/09/2023

*/

//======================= model ============================
/*
"SHE" == 1D Schrödinger-Helmholtz equation
"NLSE" == 1D nonlinear Schrödinger equation  */
const std::string FLAG_MODEL_TYPE = "SHE";


//define pi
const double pi = 3.14159265358979323846;   


//number of grid points
const int N = 512;                      

//length of periodic box
const double Lx = 2.0 * pi;                  

/* equation parameters */

//linear parameter
const double c = -0.5;

//nonlinear parameter
const double g = -1.0;

//chemical potential
const double mu = 0.0;

//SHE parameter
const double beta = 1.0;

//grid spacing
const double dx = Lx / double(N);         

//define the size of the anti-aliasing array for 3/2 - rule 
const int M = 3*N/2;                     
    
//openmp threads (not in use)
const int num_threads = 1;

//=============== timestepping ============================

/*
set timestepping method
"ETDRK4" == 4th order exponential time differencing Runge-Kutta-B 
"RK2" == 2nd order Runge-Kutta with integrating factors      */
const std::string FLAG_TIMESTEP_METHOD = "ETDRK4";

//time step
const double dt = 5.e-5;               

//total time steps
const int total_steps = 999999999;       //number of total time steps

//time at which you would like state and diagnostics to be recorded
const double output_time = 10.0;               

//set this flag to be true to record diagnostics (code will run faster if this == "false")
const bool FLAG_OUTPUT_DIAGONOSTICS = "true";

/*set initial condition for simulation
"ZERO" == generates zero state from run
"FILE" == uses initial condition from last file given curframe.dat
*/
const std::string FLAG_INITIAL_CONDITION = "ZERO";

//=========== dissipation =============================
//turn hyperviscosity on/off
const bool FLAG_HYPER_DISSIPATION = true;

//hyperviscosity coefficient
const double nu = 1.e-36;//1.e-32;	

//positive power of the laplacian (int)
const int nupower = 8;						

//turn hypoviscosity on/off
const bool FLAG_HYPO_DISSIPATION = true;

//hypoviscosity coefficient
const double alpha = 1.e3;	

//negative power of the laplacian (int) (friction = 0)				
const int alphapower = -2;				


//=================== forcing  ===============================
/*
 "off" = no forcing
"annulus" = additive forcing in annulus in fourier space:      kf-dk <= kf <= kf+dk with constant force_amplitude
"gaussian" = additive forcing using Gaussian distribution:    force_amplitude * exp (-0.5 * [(k-kf)/sigmaf)]^2 ) 
"exponential" = additive forcing using Gaussian distribution:    force_amplitude * pow(k/kf,force_power)*exp(-pow(k/kf,force_power))

*/
const std::string FLAG_FORCING_TYPE = "annulus"; //Flag to turn on additive forcing

//FLAG for having forcing amplitude rescaled to give unit energy density
const bool FLAG_FORCING_RESCALE = false; // Flag for having forcing amplitude rescaled to give unit energy density

//forcing amplitude 
const double force_amplitude = 1.e0;//0.001;	
//forcing wavenumber
const double kf = 64;	
//width of forcing annulus (in terms of wavenumber)					
const double dk = 2.0;							
//forcing standard devation for "gaussian" profile
const double sigmaf = 1.0;					
// forcing power for "exponential" forcing
const double force_power = 4.0;

//const double mean_waveaction = 1.0; //parameter for rescaling forcing wih large-scale dissipation

