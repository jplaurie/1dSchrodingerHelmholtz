#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>

#include <fftw3.h>
#include <fstream>
#include "const.h"
#include "functions.h"
using namespace std;
using namespace arma;

static double k,waveaction_inj, energy_inj;
static int forcing_modes;

void initialiseForcing(rowvec & famp){
    
if(FLAG_FORCING_TYPE == "annulus"){
	
	famp.zeros();
	waveaction_inj = 0.0;
	energy_inj = 0.0;

	for(int i = 0; i < N/2; i++){
        
		//define k
		k = abs(2.0*pi * double(-i-1)/ Lx);
    
		//define radius
		if( abs(k-kf) < dk){
			
			//set amplitude
			famp(N-i-1) = force_amplitude;
			
			//determine wave action and energy injected
			waveaction_inj += pow(famp(N-i-1),2.0);
			energy_inj += pow(famp(N-i-1),2.0)*( -c * k * k);

			//determine the number of forcing modes
			forcing_modes += 1;
		}
		
		if((i==0)) continue;                

		//define k        
		k =   2.0*pi * double(i)/ Lx;
        //define radius    
		if( abs(k-kf) <= dk){
			
			//set amplitude
			famp(i) = force_amplitude;

			//determine wave action and energy injected
			waveaction_inj += pow(famp(i),2.0);
			energy_inj += pow(famp(i),2.0)*(-c * k * k);
			//determine the number of forcing modes
			forcing_modes += 1;
		}

	

	}

}
else if(FLAG_FORCING_TYPE == "gaussian"){
	
	famp.zeros();
	waveaction_inj = 0.0;
	energy_inj = 0.0;

	for(int i = 0; i < N/2; i++){
        
		//define k
		k =  2.0*pi * double(i)/ Lx;
        //set forcing amplitude  
		famp(i) = force_amplitude * exp(-0.5 *pow( (k-kf)/sigmaf ,2.0)  );

		//determine wave action and energy injected
		waveaction_inj += pow(famp(i),2.0);
		energy_inj += pow(famp(i),2.0)*(-c * k * k);
		//determine the number of forcing modes
		forcing_modes += 1;
          
		//define k  
		k = abs(2.0*pi * double(-i-1)/ Lx);
           
		famp(N-i-1) = force_amplitude * exp(-0.5 * pow( (k-kf)/sigmaf,2.0) ) ;
		//determine wave action and energy injected
		waveaction_inj += pow(famp(N-i-1),2.0);
		energy_inj += pow(famp(N-i-1),2.0)*(-c*k*k);
		//determine the number of forcing modes
		forcing_modes += 1;
          		
	}

}
else if(FLAG_FORCING_TYPE == "exponetial"){
	
	famp.zeros();
	waveaction_inj = 0.0;
	energy_inj = 0.0;
	
	for(int i = 0; i < N/2; i++){
        
		//define k
		k =  2.0*pi * double(i)/ Lx;
        //set forcing amplitude  
		famp(i) = force_amplitude * pow(k/kf,force_power)*exp(-pow(k/kf,force_power));

		//determine wave action and energy injected
		waveaction_inj += pow(famp(i),2.0);
		energy_inj += pow(famp(i),2.0)*(-c*k*k);
		//determine the number of forcing modes
		forcing_modes += 1;
          
		//define k  
		k = abs(2.0*pi * double(-i-1)/ Lx);
           
		famp(N-i-1) = force_amplitude * pow(k/kf,force_power)*exp(-pow(k/kf,force_power));
		
		//determine wave action and energy injected
		waveaction_inj += pow(famp(N-i-1),2.0);
		energy_inj += pow(famp(N-i-1),2.0)*(-c*k*k);
		//determine the number of forcing modes
		forcing_modes += 1;
          		
	}
}
else{
	cout << "FLAG_FORCING_TYPE not defined" << endl;
	exit(1);
}

    
if(FLAG_FORCING_RESCALE == true){
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



