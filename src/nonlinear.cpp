#include <armadillo>
#include <fftw3.h>
#include <complex>
#include <cmath>
#include "const.h"

using namespace std;
using namespace arma;

static cx_rowvec Npsi_temp_M(M,fill::zeros), psi_hat_N(N,fill::zeros);
static cx_rowvec k2(M,fill::zeros);

//function conditioners for embedding data in larger arrays anti-aliasing
void embed_N_M( cx_rowvec, cx_rowvec &);
void embed_M_N( cx_rowvec, cx_rowvec &);
void dealias(cx_rowvec &);
void define_k2(cx_rowvec &);


cx_rowvec nonlinear( cx_rowvec A , double da, cx_rowvec L, cx_rowvec & psi_M,cx_rowvec & psi_hat_M, fftw_plan FFT, fftw_plan IFFT){ 
/*
	Subroutine to compute the nonlinear term in the 1S SHE.NLS equation.  Use a 3/2 dealias rule but applied twice to preserve momentum conservation

*/
//set input psi_hat = A to array got FFT
psi_hat_N.zeros();
psi_hat_N = A;

//define k2 array
define_k2(k2);

//if using RK2 apply integrating factor
if(FLAG_TIMESTEP_METHOD == "RK2"){
	psi_hat_N %= exp(da*L);
}
    
//set arrays to zero
Npsi_temp_M.zeros();
psi_hat_M.zeros();
psi_M.zeros();
   
//embeds data in larger array
embed_N_M(psi_hat_N,psi_hat_M);             

//inverse fft (size M)
fftw_execute(IFFT);                 

//temporally store psi_M
Npsi_temp_M = psi_M;	

//multiply by conjugate
psi_M %= conj(psi_M);
 	
//forward fft (size M)
fftw_execute(FFT);

//normalize
psi_hat_M /= double(M);            

//compute nonlinear term based on model
if(FLAG_MODEL_TYPE == "SHE"){  
    psi_hat_M %= g*pow(1.0+(beta_param*k2), -1.0);
}
else if(FLAG_MODEL_TYPE == "NLSE"){  
    psi_hat_M *= g;
}
else{
    cout << "FLAG_MODEL_TYPE invalid" << endl;
    exit(1);
}

//dealias using 3/2 rule
dealias(psi_hat_M);

//inverse fft (size M)
fftw_execute(IFFT);

//muliply by stored |PSI|^2
psi_M %= Npsi_temp_M;
 
//forward fft (size M)
fftw_execute(FFT);                  

//normalize
psi_hat_M /= double(M);    

//put data back into smaller array size N
embed_M_N(psi_hat_M,psi_hat_N);            

//divide by i (from time derivative)
psi_hat_N /= complex<double>(0.0,1.0) ;

//apply integrating factor if using RK2
if(FLAG_TIMESTEP_METHOD == "RK2"){
	psi_hat_N %= exp(-da*L);
}

//return nonlinear term
return psi_hat_N;                   
}



//================================================================================================================
void embed_N_M( cx_rowvec A, cx_rowvec & B){
/* embeds data of size N into size M array */    

    B.zeros();
    
    for(int i =0; i < N/2 ; i++){    
        B(i) = A(i);
        B(M-i-1)= A(N-i-1);       
    }
  
    return;
}


void embed_M_N( cx_rowvec B, cx_rowvec & A){
/* embeds data of size M into size N array */ 
    
    A.zeros();
    
    for(int i =0; i < N/2 ; i++){
        A(i) = B(i);
        A(N-i-1)= B(M-i-1);
    }   
    return;
}
    
 
void dealias(cx_rowvec & A){
/*dealiases size M array using 3/2 rule*/
    for(int i = N/2; i < M - (N/2) ; i++){
        A(i) = complex<double>(0.0,0.0);
    }

    return;
}
       
void define_k2(cx_rowvec & k2_temp){
     /* defines k2 array */
    for(int i = 0; i < M/2 ; i++){
        k2_temp(i) = complex<double>(pow(2.0 * pi * double(i) / Lx,2.0),0.0);
        k2_temp(M-i-1) = complex<double>(pow(2.0 * pi * double(-i-1) / Lx,2.0),0.0);
    }

return;
}
    
    
