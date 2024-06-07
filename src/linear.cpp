//this creates the linear operator in fourier-space


#include <armadillo>
#include <fftw3.h>
#include "const.h"
#include "functions.h"
#include <cmath>
#include <complex>
using namespace std;
using namespace arma;


void linearOperator(cx_rowvec & L){ 
    /*
     subroutine to define the linear operator of the laplacian and dissipation
    */


    //define linear derivative and chemical potential
	for(int i = 0; i < N/2; i++){
		
        L(i) = complex<double>( - c  * pow(2.0 * pi * double(i) / Lx,2.0)    + mu , 0.0 );
		L(N-i-1) = complex<double>(- c  * pow(2.0 * pi * double(-i-1) / Lx,2.0) + mu , 0.0 );
		
	}

    //divide by i (from time derivative)
    L /= complex<double>(0.0,1.0);


    //add contribution from hyperviscosity
    if(FLAG_HYPER_DISSIPATION == true){

        for(int i = 0; i < N/2; i++){
            L(i) -= complex<double>( nu * pow(   pow(  2.0 * pi * double(i) /Lx ,2.0)  , nupower ),0.0);
            L(N-i-1) -= complex<double>( nu * pow(pow((2.0 * pi * double(-i-1) / Lx),2.0) , nupower ),0.0);
        }
    }

    //add contribution from hypoviscosity
    if(FLAG_HYPO_DISSIPATION == true){

        for(int i = 0; i < N/2; i++){
            L(i) -=  complex<double>( alpha * pow(pow(2.0 * pi * double(i) / Lx,2.0), alphapower),  0.0);
            L(N-i-1) -=  complex<double>( alpha * pow(pow((2.0 * pi * double(-i-1) / Lx),2.0), alphapower),0.0);
    
        }
        //removes zeroth mode contirbution in case alphapower < 0
        if(alphapower < 0){
            L(0) = complex<double>(0.0,0.0);
        } 

    }
return;
}



void setupETDRK(cx_rowvec L_E,cx_rowvec & Q1,cx_rowvec & Q2,cx_rowvec & Q3,cx_rowvec & Q4,cx_rowvec & Q5, cx_rowvec & F1, cx_rowvec & F2, cx_rowvec & F3, cx_rowvec & E1, cx_rowvec & E2, double h){

    /*
    Routine to set up commonly used operators for ETDRK4 timestepping routine. This is the code for the original ETDRK4 from Cox and Matthews and also ETDRK-B (same)
    */

    cx_rowvec LR(N,fill::zeros);
    
    //number of complex contour points
    int M = 32;									

    E1 = exp(h*L_E);
    E2 = exp(0.5*h*L_E);

    for(int i  = 0; i < M; i++){

        LR = h*L_E + exp(complex<double>(0.0,1.0)*2.0*pi*(double(i)-0.5)/ double(M));
    
        Q1 += (exp(0.5*LR) - 1.0) / LR;
        Q2 += (exp(0.5*LR)%(LR-4.0) + LR + 4.0) / pow(LR,2.0); 
        Q3 += 2.0*(2.0*exp(0.5*LR) - LR - 2.0) / pow(LR,2.0); 
        Q4 += (exp(LR)%(LR-2.0) + LR + 2.0) / pow(LR,2.0); 
        Q5 += 2.0*(exp(LR) - LR - 1.0) / pow(LR,2.0); 
        F1 += (-4.0 - LR + exp(LR)%(4.0-3.*LR+pow(LR,2.0))) / pow(LR,3.0);
        F2 += (2.0 + LR + exp(LR)%(-2.0 + LR)) / pow(LR,3.0);
        F3 += (-4.0 - 3.0*LR - pow(LR,2.0) + exp(LR)%(4.0-LR)) / pow(LR,3.0);
    }

    Q1 *= h / double(M);
    Q2 *= h / double(M);
    Q3 *= h / double(M);
    Q4 *= h / double(M);
    Q5 *= h / double(M);
    F1 *= h / double(M);
    F2 *= h / double(M);
    F3 *= h / double(M);

    return;
}
