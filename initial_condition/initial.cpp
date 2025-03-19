
//standalone file that creates an initial condition

#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <chrono>

using namespace std;
using namespace arma;

const int N = 4096;
const double pi = 3.14159265358979323846;
const double L = 2.0*pi;
const double dx= L/double(N);

const bool FLAG_INITIAL_FOURIER = true;
const double kf = 16.0;
const double sigma_k = 4.0;

const double amp = 10.0;

int main(){
	
    cx_rowvec psi(N,fill::zeros);
    cx_rowvec psi_hat(N,fill::zeros);

if(FLAG_INITIAL_FOURIER == true){
double phase = 0.0;

    fftw_plan IFFT;
    IFFT = fftw_plan_dft_1d(N, (fftw_complex*) psi_hat.memptr(), (fftw_complex*) psi.memptr(), FFTW_BACKWARD, FFTW_PATIENT);

    //set up random forcing
    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();  
    mt19937_64 generator(seed1);
    uniform_real_distribution<double> distribution(0.0,2.0*pi);
    
    double k=0.0;
    
    for(int i =0; i< N/2; i++){

        k =  2.0*pi * double(i)/ L;
        phase =  distribution(generator);
        psi_hat(i) = amp*exp(-0.5*pow( (k - kf)/sigma_k ,2.0))*complex<double>(cos(phase),sin(phase));

        
        k =abs(2.0*pi * double(-i-1)/ L);
        phase =  distribution(generator);
        psi_hat(N-i-1) = amp*exp(-0.5*pow( (k - kf)/sigma_k ,2.0))*complex<double>(cos(phase),sin(phase));

    }


    fftw_execute(IFFT);    


}
else{

    double A= 1.0;

    for(int i=0; i < N; i++){ 
       // psi(i) = A*pow(cosh( dx *(i - (N/2))),-1.0)*complex<double>(0.7071,0.7071);
       // psi(i) = A*cos(10 * i*dx) + A*sin(11.*i*dx) + A*cos(12.*i*dx);
	 psi(i) = complex<double>(A,0.0);    
}

    
}

ofstream fout("./psi.000000");              //open up file
    fout.precision(12);
    for(int i=0; i < N; i++){
        
        fout << i*dx << " " << real(psi(i)) << " " << imag(psi(i)) << endl;
       
    }

    fout.close();

    return 0;

}
