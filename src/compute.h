void computeEnergy(cx_rowvec , cx_rowvec , double & , double & , double & ,cx_rowvec & ,cx_rowvec & ,fftw_plan , fftw_plan );
void computeWaveaction(cx_rowvec , double & );
void computeDissipationRate( cx_rowvec , double & , double & , double & , double & );
void computeSpectrum( cx_rowvec , rowvec & );
void computeFlux(cx_rowvec , cx_rowvec L,cx_rowvec & ,cx_rowvec &, rowvec &  ,rowvec & ,fftw_plan , fftw_plan );