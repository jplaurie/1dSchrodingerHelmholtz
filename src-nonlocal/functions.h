
//this includes all the preconditoning plans for the funtions


#include <armadillo>
#include <complex>
#include <fftw3.h>

using namespace std;
using namespace arma;

void linear(cx_rowvec &);

cx_rowvec nonlinear(cx_rowvec , double, cx_rowvec ,cx_rowvec &,cx_rowvec &,fftw_plan,fftw_plan);


void read(cx_rowvec &,cx_rowvec &, double& , int&, fftw_plan);
void energy(double, int,cx_rowvec , cx_rowvec, double &,cx_rowvec & ,cx_rowvec &,fftw_plan, fftw_plan);

void intensity(double, int, cx_rowvec );

void spectrum( cx_rowvec, rowvec &, rowvec & ,int, int );
void flux( cx_rowvec,cx_rowvec, rowvec &, rowvec & ,int, int);
void forcing(cx_rowvec,rowvec);

void initial_forcing1(rowvec & , double, double);
void initial_forcing2(rowvec & , double, double);

void record_parameters();

void setup_ETDRK(cx_rowvec,cx_rowvec &,cx_rowvec &,cx_rowvec &,cx_rowvec &, cx_rowvec &, cx_rowvec &,cx_rowvec &, cx_rowvec &, cx_rowvec &, cx_rowvec &, double);