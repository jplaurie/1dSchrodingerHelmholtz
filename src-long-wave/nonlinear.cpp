#include <armadillo>
#include <fftw3.h>
#include "const.h"
#include "functions.h"
#include <complex>
#include <cmath>

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
	Subroutine to compute the nonlinear cubic term in the full 1D OWT equation.  Use a 3/2 dealias rule but applied twice to preserve momentum conservation

*/

psi_hat_N.zeros();
psi_hat_N = A;

define_k2(k2);


if(FLAG_TIMESTEP_ETDRK == false){
	psi_hat_N %= exp(da*L);
}
    
Npsi_temp_M.zeros();
psi_hat_M.zeros();
psi_M.zeros();
   
embed_N_M(psi_hat_N,psi_hat_M);             //embeds data in larger array
  
fftw_execute(IFFT);                  //inverse fft size M

Npsi_temp_M = psi_M;	//temp stores psi_M

psi_M %= conj(psi_M);
 	
fftw_execute(FFT);

psi_hat_M /= double(M);            //normalize

psi_hat_M %= -0.5*pow(1.0+(k2/g), -1.0);

dealias(psi_hat_M);

fftw_execute(IFFT);

psi_M %= Npsi_temp_M;
 
fftw_execute(FFT);                  //forward fft size M

psi_hat_M /= double(M);            //normalize

embed_M_N(psi_hat_M,psi_hat_N);             //puts data back into smaller array size N

psi_hat_N /= complex<double>(0.0,1.0) ;

if(FLAG_TIMESTEP_ETDRK == false){
	psi_hat_N %= exp(-da*L);
}

return psi_hat_N;                   //returns nonlinear term
}



//================================================================================================================
void embed_N_M( cx_rowvec A, cx_rowvec & B){
    
    B.zeros();
    
        for(int i =0; i < N/2 ; i++){    
            B(i) = A(i);
            B(M-i-1)= A(N-i-1);       
        }
  
    return;
}


void embed_M_N( cx_rowvec B, cx_rowvec & A){
    
A.zeros();
    
    for(int i =0; i < N/2 ; i++){
        A(i) = B(i);
        A(N-i-1)= B(M-i-1);
    }

    
return;
}
    
 
void dealias(cx_rowvec & A){
     
    for(int i = N/2; i < M - (N/2) ; i++){
        A(i) = complex<double>(0.0,0.0);
    }

return;
}
       
void define_k2(cx_rowvec & k2_temp){
     
    for(int i = 0; i < M/2 ; i++){
        k2_temp(i) = complex<double>(pow(2.0 * pi * double(i) / Lx,2.0),0.0);
        k2_temp(M-i-1) = complex<double>(pow(2.0 * pi * double(-i-1) / Lx,2.0),0.0);
    }

return;
}
    
    
