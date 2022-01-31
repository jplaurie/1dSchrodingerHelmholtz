
#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include "const.h"
#include "functions.h"
#include <fftw3.h>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace arma;

void intensity(double runtime, int out_count, cx_rowvec psi){
    
ostringstream out_I;
out_I << "./data/psi." << setw(6) << setfill('0') << out_count << ends;			//creates file name for outputting data at time slice
string filename = out_I.str();
ofstream fout_I(filename.c_str());
fout_I << scientific;
fout_I.precision(12);
    
for(int i=0; i < N; i++){
		fout_I << i*dx << " " << real(psi(i)) << " " << imag(psi(i)) << endl;
	}


fout_I.close();

return;
    
}