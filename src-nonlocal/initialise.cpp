#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include "const.h"
#include "functions.h"
#include <fftw3.h>
#include <fstream>

using namespace std;
using namespace arma;

static double k,waveaction_inj, energy_inj;
static int forcing_modes;

void initial_forcing1(rowvec & famp, double kf, double dk){
    
	for(int i =0; i < N/2; i++){
        
		
		k = abs(2.0*pi * double(-i-1)/ Lx);
               	 
		if( abs(k-kf) < dk){
			
			famp(N-i-1) = Amp;
			waveaction_inj += pow(famp(N-i-1),2.0);
			energy_inj += pow(famp(N-i-1),2.0)*(k*k);
			forcing_modes += 1;
		}
		
		if((i==0)) continue;                
                
		k =   2.0*pi * double(i)/ Lx;
            
		if( abs(k-kf) < dk){
			 
			famp(i) = Amp;
			waveaction_inj += pow(famp(i),2.0);
			energy_inj += pow(famp(i),2.0)*(k*k);
			forcing_modes += 1;
		}

	

	}

    
if(FLAG_FORCE_AMP_RESCALE == true){
	famp *= sqrt(2.0*alpha/waveaction_inj);
	energy_inj *= (2.0*alpha/waveaction_inj);
	waveaction_inj = 2.0*alpha;
}
    
cout << "number of forcing modes = " << forcing_modes << endl;
cout << "energy injection rate = " << energy_inj << endl;
cout << "waveaction injection rate = " << waveaction_inj << endl;
    	
ofstream fout_force_data("./output/force_values.txt");
fout_force_data.precision(12);

fout_force_data << "number of forcing modes = " << forcing_modes << endl;
fout_force_data << "energy injection rate = " << energy_inj << endl;
fout_force_data << "waveaction injection rate = " << waveaction_inj << endl;
fout_force_data << "forcing mode kf = " << kf << endl;
fout_force_data << "forcing width dk = " << dk << endl;

fout_force_data.close();
    
ofstream fout_Amp("./output/force_spectrum.dat");
fout_Amp << scientific;
fout_Amp.precision(12);
 
for(int i =N/2; i< N; i++){
       
	k = 2.0*pi * double(i-N)/ Lx;
	fout_Amp << k << " " << famp(i) << endl;
}
	
for(int i =0; i< N/2; i++){
    		
	k = 2.0*pi * double(i)/ Lx;
	fout_Amp << k << " " << famp(i) << endl;     
}
	

fout_Amp.close();

return;
}



void initial_forcing2(mat & famp, double kf, double dk){
    
    
forcing_modes =0;
famp=0.0;
energy_inj=0.0;
waveaction_inj=0.0;
    
	for(int i =0; i< N/2; i++){

		k =  2.0*pi * double(i)/ Lx;
            
		famp(i) = pow(k/kf,force_power)*exp(-pow(k/kf,force_power));
		waveaction_inj += pow(famp(i),2.0);
		energy_inj += pow(famp(i),2.0)*(k*k);
		forcing_modes += 1;
            
		k =abs(2.0*pi * double(-i-1)/ Lx);
           
		famp(N-i-1) = pow(k/kf,force_power)*exp(-pow(k/kf,force_power));
		waveaction_inj += pow(famp(N-i-1),2.0);
		energy_inj += pow(famp(N-i-1),2.0)*(k*k);
		forcing_modes += 1;
          		
	}

    
if(FLAG_FORCE_AMP_RESCALE == true){
	famp *= sqrt(2.0*alpha/waveaction_inj);
	energy_inj *= (2.0*alpha/waveaction_inj);
	waveaction_inj = 2.0*alpha;
}

cout << "number of forcing modes = " << forcing_modes << endl;
cout << "energy injection rate = " << energy_inj << endl;
cout << "waveaction injection rate = " << waveaction_inj << endl;
        
ofstream fout_force_data("./output/force_values.txt");
fout_force_data.precision(12);

fout_force_data << "number of forcing modes = " << forcing_modes << endl;
fout_force_data << "energy injection rate = " << energy_inj << endl;
fout_force_data << "waveaction injection rate = " << waveaction_inj << endl;
fout_force_data << "forcing mode kf = " << kf << endl;
fout_force_data << "forcing width dk = " << dk << endl;
fout_force_data.close();

ofstream fout_Amp("./output/force_spectrum.dat");
fout_Amp << scientific;
fout_Amp.precision(12);

for(int i =N/2; i< N; i++){
       
	k = 2.0*pi * double(i-N)/ Lx;


		fout_Amp << k << " " << famp(i) << endl;
	}

for(int i =0; i< N/2; i++){
    		
	k = 2.0*pi * double(i)/ Lx;
       
		fout_Amp << k << " " << famp(i) << endl;     
	}
        
    
fout_Amp.close();
    
return;
}



