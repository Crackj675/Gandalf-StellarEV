//=================================================================================================
//  BonnorEbertSophereIc.cpp
//  Class for generating initial conditions for ...
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=================================================================================================


#include <fstream>
#include <sstream>
#include "Precision.h"
#include "Debug.h"
#include "Ic.h"
using namespace std;

//=================================================================================================
//  Ic::Bonnor-Eber Sphere
/// Set-up Bonnor-Eber Sphere (1956) initial conditions
//=================================================================================================
template <int ndim>
BonnorEbertSphereIc<ndim>::BonnorEbertSphereIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
}


template <int ndim>
void BonnorEbertSphereIc<ndim>::Generate(void)
{
  string ic = simparams->stringparams["ic"];
  if (ic == "bes") {
		int i;                               // Particle counter
		int j;
		int k;                               // Dimension counter
		int Nsphere;                         // Actual number of particles in sphere
		int NHotGas;
		int Nshell;
		FLOAT gpecloud;                      // Total grav. potential energy of entire cloud
		FLOAT keturb;                        // Total turbulent kinetic energy of entire cloud
		FLOAT mp;                            // Mass of one particle
		FLOAT rcentre[ndim];                 // Position of sphere centre
		FLOAT rho;                           // Fluid density
		FLOAT xmin;                          // Minimum coordinate value
		FLOAT vfactor;                       // Velocity scaling factor (to scale correct alpha_turb)
		FLOAT *r;                            // Positions of all particles
		FLOAT *r_shell;                            // Positions of all particles
		FLOAT *v;                            // Velocities of all particles
		FLOAT *rHotGas;                            // Velocities of all particles
		FLOAT dxgrid;                        // Grid spacing
		FLOAT rmax[ndim];                    // Maximum size of bounding box
		FLOAT rmin[ndim];                    // Minimum size of bounding box
		FLOAT Cs;							   // Sound speed
		FLOAT r0;
		FLOAT r_bar;
		FLOAT M_shell;
		FLOAT mcore;
		FLOAT Etherm;
		DOUBLE *vfield;                      // Table with turbulent velocity field from


		// Create local copies of initial conditions parameters
		int field_type   		= simparams->intparams["field_type"];
		int gridsize     		= simparams->intparams["gridsize"];
		int Npart        		= simparams->intparams["Nhydro"];
		FLOAT alpha_turb 		= simparams->floatparams["alpha_turb"];
		FLOAT gammaone   		= simparams->floatparams["gamma_eos"] - 1.0;
		FLOAT mcloud     		= simparams->floatparams["mcloud"];
		FLOAT mu_bar     		= simparams->floatparams["mu_bar"];
		FLOAT power_turb 		= simparams->floatparams["power_turb"];
		FLOAT radius    		= simparams->floatparams["radius"];
		FLOAT temp0      		= simparams->floatparams["temp0"];
		FLOAT rho_ambient   = simparams->floatparams["rho_ambient"];
		FLOAT T_ambient	   	= simparams->floatparams["T_ambient"];
		FLOAT R_ambient	   	= simparams->floatparams["R_ambient"];
		FLOAT alpha_vir     = simparams->floatparams["alpha_vir"];
		FLOAT fkmax 			  = simparams->floatparams["fkmax"];
		string line;

		// Convert any parameters to code units
		mcloud 		/= simunits.m.outscale;
		radius 		/= simunits.r.outscale;
		temp0  		/= simunits.temp.outscale;
		rho_ambient	/= simunits.rho.outscale;
		T_ambient	    /= simunits.temp.outscale;
		R_ambient	    /= simunits.r.outscale;



	#if !defined(FFTW_TURBULENCE)
		string message = "FFTW turbulence flag not set";
		ExceptionHandler::getIstance().raise(message);
	#endif

		debug2("[Ic::Bonnor-Eber Sphere]");

		cout << "Creating Bonnor-Ebert Sphere" << endl;



		// Read in positions of Particle from file

		//fileName = param->filePath + sOutFileID;
		ifstream pFile ("BESpos.dat");
		if (pFile.is_open()){
			getline(pFile, line);
			stringstream ss(line);
			ss >> NHotGas;
			getline(pFile, line);
			stringstream ss1(line);
			ss1 >> Npart;
			getline(pFile, line);
			stringstream ss2(line);
			ss2 >> r0;
			r0 /= simunits.r.outscale;
			radius = 6.5 *r0;
			cout << "Npart = " << Npart << " NHotGas = " << NHotGas << " r0 = " << r0 << endl;
			//r0 = (r0/(100*3.0857*pow(10,16))) / simunits.r.outscale;

			rHotGas = new FLOAT[ndim*NHotGas];
			r = new FLOAT[ndim*Npart];
			v = new FLOAT[ndim*Npart];

			int i;
			int j;
			i = 0;
			j = 0;

			// read in psositions
			while(getline(pFile, line) and i < Npart){
				stringstream ss(line);
				ss >> r[ndim*i] >> r[ndim*i+1] >> r[ndim*i+2] >> mp;
				for (k=0; k<ndim; k++) r[ndim*i + k] *= r0;
					i++;
				}
				while(getline(pFile, line) and j < NHotGas+1){
					stringstream ss(line);
					ss >> rHotGas[ndim*j] >> rHotGas[ndim*j+1] >> rHotGas[ndim*j+2] >> mp;

					for (k=0; k<ndim; k++) rHotGas[ndim*j + k] *= r0;
		          j++;
					}
		     pFile.close();
		  }
		  else{
		    cout << "Unable to open file";
		  }

		// make BES super critical
		mp *= 2.0;

		for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;


		M_shell = (4.0*onethird*pi*pow(R_ambient,3) - 4.0*onethird*pi*pow(6.5 * r0,3)) * rho_ambient;
		Nshell = int (M_shell / mp);
		r_shell = new FLOAT[ndim*Nshell];
		cout << "Ncore: " << Npart << "   NhotGas:  " << NHotGas << "   Nshell:  " << Nshell << "   Mshell: " << M_shell << endl; 
		Ic<ndim>::SphericalShell(Nshell, rcentre, radius, R_ambient, r_shell);


		// Allocate local and main particle memory
		hydro->Nhydro = Npart+NHotGas-1 + Nshell;
		sim->AllocateParticleMemory();
		mcloud = mp * (Npart + NHotGas + Nshell);
		mcore =  mp * (Npart);
		cout << "Mcloud = " << mcloud << " Mcore = " << mcore << endl;
		rho = (FLOAT) 3.0*mcore / ((FLOAT) 4.0*pi*pow(radius,3));

		j = 0;
		// Record particle properties in main memory
		for (i=0; i<hydro->Nhydro; i++) {
			Particle<ndim>& part = hydro->GetParticlePointer(i);
			if (i < Npart){
				for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
				for (k=0; k<ndim; k++) part.v[k] = 0.0;
				part.m = mp;
				part.u = temp0/gammaone/mu_bar;
				part.temp_ambient = temp0;		
				part.outflowFlag = (int) 0;
		   	part.iorig = i;	
			}
			if (i >= Npart and i < Npart+NHotGas){
		 		r_bar = 0;
				for (k=0; k<ndim; k++) r_bar += pow(rHotGas[ndim*j+k],2);
				r_bar = sqrt(r_bar);
				if (r_bar >= 6.5*r0){
					for (k=0; k<ndim; k++) part.r[k] = rHotGas[ndim*j+k];
					for (k=0; k<ndim; k++) part.v[k] = 0.0;
					part.m = mp;
					part.u = (temp0)/gammaone/mu_bar;
					part.temp_ambient = temp0;		
					part.outflowFlag = (int) 2;
					j++;
				}
				else{
					for (k=0; k<ndim; k++) part.r[k] = rHotGas[ndim*j + k];
					for (k=0; k<ndim; k++) part.v[k] = 0.0;
					part.m = mp;
					part.u = temp0/gammaone/mu_bar;
					part.temp_ambient = temp0;		
					part.outflowFlag = (int) 0;
					j++;
				}
			}
			if(i >= Npart+NHotGas){
		  	for (k=0; k<ndim; k++) part.r[k] = r_shell[ndim*(i-(Npart+NHotGas)) + k];
		    for (k=0; k<ndim; k++) part.a[k] = (FLOAT) 0.0;	
				part.u = T_ambient/gammaone/mu_bar;
				part.temp_ambient = temp0;		
				part.m = mp;
				part.outflowFlag = (int) 2;
			}
		}
		
		sim->initial_h_provided = false;


		cout << "kmax "<< fkmax << endl;
		// Generate turbulent velocity field for given power spectrum slope



		vfield = new DOUBLE[ndim*gridsize*gridsize*gridsize];
		dxgrid = 2*radius/(FLOAT) (gridsize - 1);
		xmin = -radius;

		// Generate gridded velocity field
		Ic<ndim>::GenerateTurbulentVelocityField(field_type, gridsize, power_turb, vfield);
	
		// Now interpolate generated field onto particle positions
		Ic<ndim>::InterpolateVelocityField(Npart, gridsize, xmin, dxgrid, r, vfield, v);

		// Finally, copy velocities to main SPH particle array
		for (i=0; i<Npart+NHotGas-1; i++) {
			Particle<ndim>& part = hydro->GetParticlePointer(i);
			for (k=0; k<ndim; k++) part.v[k] = v[ndim*i + k];
		}

		// Change to COM frame of reference
		sim->SetComFrame();

		// Calculate sound speed
		Cs = pow( (gammaone+1)*(temp0/mu_bar)*(k_boltzmann/m_hydrogen) / (simunits.v.outscale*simunits.v.outSI*simunits.v.outscale*simunits.v.outSI/(simunits.temp.outscale*simunits.temp.outSI)),0.5);
		Etherm = 0.5*Cs*Cs*mcore;
		gpecloud = (FLOAT) ((FLOAT) 3 / (FLOAT) 5) *mcore*mcore/(radius);



		// Calculate total kinetic energy of turbulent velocity field
		keturb = (FLOAT) 0.0;
		for (i=0; i<Npart+NHotGas-1; i++) {
		  Particle<ndim>& part = hydro->GetParticlePointer(i);
			for (k=0; k<ndim; k++) keturb += part.v[k]*part.v[k];
		}

		// Resacel velocities to match alpha_vir
		keturb *= (FLOAT) mp*0.5;
		cout << "ketturb: " << keturb << "  therm: " << Etherm << " Egrav: " << gpecloud << " alpha: " << alpha_vir << " mcloud: " << mcloud <<  endl;
		vfactor = sqrt((gpecloud*alpha_vir - 2*Etherm)/(2*keturb));


		for (i=0; i<Npart+NHotGas-1; i++) {
		  Particle<ndim>& part = hydro->GetParticlePointer(i);
		  for (k=0; k<ndim; k++) part.v[k] *= vfactor;
		}

		// Calculate total kinetic energy of turbulent velocity field
		keturb = (FLOAT) 0.0;
		for (i=0; i<Npart+NHotGas-1; i++) {
		  Particle<ndim>& part = hydro->GetParticlePointer(i);
			for (k=0; k<ndim; k++) keturb += part.v[k]*part.v[k];
		}
		keturb *= (FLOAT) mp*0.5;

		cout << "ketturb: " << keturb << "  therm: " << Etherm << " Egrav: " << gpecloud << " vfactor: " << vfactor <<  endl;
		cout << "alpha_vir : " << (2*Etherm + 2*keturb)/gpecloud  <<" --- " << alpha_vir <<  "   Epot/Ekin : " << gpecloud/(keturb+Etherm) 	<< "    Machnumber : " << pow((2*keturb/mcore), 0.5)/ Cs << "    Sound speed : " << Cs*simunits.v.outscale << endl;
	//(Cs*Cs + (2*keturb/mcloud))*radius/(mcloud)

		for (i=Npart+NHotGas-1; i<hydro->Nhydro; i++) {
			Particle<ndim>& part = hydro->GetParticlePointer(i);
			for (k=0; k<ndim; k++) part.v[k] = 0.0;		
		}

		delete[] v;
		delete[] r;
		delete[] rHotGas;
		delete[] r_shell;

		return;
	}

	if (ic == "besrot" || ic == "BesRot") {

		int i;                               // Particle counter
		int j;
		int k;                               // Dimension counter
		int Nsphere;                         // Actual number of particles in sphere
		int NHotGas;
		int Nshell;
		int Npart;
		FLOAT mp;                            // Mass of one particle
		FLOAT rcentre[ndim];                 // Position of sphere centre
		FLOAT rho;                           // Fluid density
		FLOAT xmin;                          // Minimum coordinate value
		FLOAT vfactor;                       // Velocity scaling factor (to scale correct alpha_turb)
		FLOAT r0;
		FLOAT *v;                            // Velocities of all particles
		FLOAT *r;                            // Positions of all particles
		FLOAT *r_shell;                            // Positions of all particles
		FLOAT *rHotGas;                            // Velocities of all particles
    FLOAT r_perp[3];                     // Perpendicular vector
    FLOAT pc = 3.08567758*pow(10,13);    // ..

		FLOAT r_bar;
		FLOAT M_shell;
		FLOAT mcore;
		FLOAT Etherm;
		DOUBLE *vfield;                      // Table with turbulent velocity field from
	

		// Create local copies of initial conditions parameters
		FLOAT gammaone   		= simparams->floatparams["gamma_eos"] - 1.0;
		FLOAT mcloud     		= simparams->floatparams["mcloud"];
		FLOAT mu_bar     		= simparams->floatparams["mu_bar"];
		FLOAT radius    		= simparams->floatparams["radius"];
		FLOAT temp0      		= simparams->floatparams["temp0"];
		FLOAT rho_ambient   = simparams->floatparams["rho_ambient"];
		FLOAT T_ambient	   	= simparams->floatparams["T_ambient"];
		FLOAT R_ambient	   	= simparams->floatparams["R_ambient"];
    FLOAT omega          = simparams->floatparams["omega"];
		string line;

		// Convert any parameters to code units
		mcloud 		/= simunits.m.outscale;
		radius 		/= simunits.r.outscale;
		temp0  		/= simunits.temp.outscale;
		rho_ambient	/= simunits.rho.outscale;
		T_ambient	    /= simunits.temp.outscale;
		R_ambient	    /= simunits.r.outscale;



	#if !defined(FFTW_TURBULENCE)
		string message = "FFTW turbulence flag not set";
		ExceptionHandler::getIstance().raise(message);
	#endif

		debug2("[Ic::Bonnor-Eber Sphere]");

		cout << "Creating Bonnor-Ebert Sphere" << endl;



		// Read in positions of Particle from file

		//fileName = param->filePath + sOutFileID;
		ifstream pFile ("BESpos.dat");
		if (pFile.is_open()){
			getline(pFile, line);
			stringstream ss(line);
			ss >> NHotGas;
			getline(pFile, line);
			stringstream ss1(line);
			ss1 >> Npart;
			getline(pFile, line);
			stringstream ss2(line);
			ss2 >> r0;
			r0 /= simunits.r.outscale;
			radius = 6.5 *r0;
			cout << "Npart = " << Npart << " NHotGas = " << NHotGas << " r0 = " << r0 << endl;
			//r0 = (r0/(100*3.0857*pow(10,16))) / simunits.r.outscale;

			rHotGas = new FLOAT[ndim*NHotGas];
			r = new FLOAT[ndim*Npart];
			v = new FLOAT[ndim*Npart];

			int i;
			int j;
			i = 0;
			j = 0;

			// read in psositions
			while(getline(pFile, line) and i < Npart){
				stringstream ss(line);
				ss >> r[ndim*i] >> r[ndim*i+1] >> r[ndim*i+2] >> mp;
				for (k=0; k<ndim; k++) r[ndim*i + k] *= r0;
					i++;
				}
				while(getline(pFile, line) and j < NHotGas+1){
					stringstream ss(line);
					ss >> rHotGas[ndim*j] >> rHotGas[ndim*j+1] >> rHotGas[ndim*j+2] >> mp;

					for (k=0; k<ndim; k++) rHotGas[ndim*j + k] *= r0;
		          j++;
					}
		     pFile.close();
		  }
		  else{
		    cout << "Unable to open file";
		  }

		// make BES super critical
		mp *= 2.0;

		for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;


		M_shell = (4.0*onethird*pi*pow(R_ambient,3) - 4.0*onethird*pi*pow(6.5 * r0,3)) * rho_ambient;
		Nshell = int (M_shell / mp);
		r_shell = new FLOAT[ndim*Nshell];
		cout << "Ncore: " << Npart << "   NhotGas:  " << NHotGas << "   Nshell:  " << Nshell << "   Mshell: " << M_shell << endl; 
		Ic<ndim>::SphericalShell(Nshell, rcentre, radius, R_ambient, r_shell);


		// Allocate local and main particle memory
		hydro->Nhydro = Npart+NHotGas-1 + Nshell;
		sim->AllocateParticleMemory();
		mcloud = mp * (Npart + NHotGas + Nshell);
		mcore =  mp * (Npart);
		cout << "Mcloud = " << mcloud << " Mcore = " << mcore << endl;
		rho = (FLOAT) 3.0*mcore / ((FLOAT) 4.0*pi*pow(radius,3));

		j = 0;
		// Record particle properties in main memory
		for (i=0; i<hydro->Nhydro; i++) {
			Particle<ndim>& part = hydro->GetParticlePointer(i);
			if (i < Npart){
				for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
				for (k=0; k<ndim; k++) part.v[k] = 0.0;
				part.m = mp;
				part.u = temp0/gammaone/mu_bar;
				part.temp_ambient = temp0;		
				part.outflowFlag = (int) 0;
		   	part.iorig = i;	
			}
			if (i >= Npart and i < Npart+NHotGas){
		 		r_bar = 0;
				for (k=0; k<ndim; k++) r_bar += pow(rHotGas[ndim*j+k],2);
				r_bar = sqrt(r_bar);
				if (r_bar >= 6.5*r0){
					for (k=0; k<ndim; k++) part.r[k] = rHotGas[ndim*j+k];
					for (k=0; k<ndim; k++) part.v[k] = 0.0;
					part.m = mp;
					part.u = (temp0)/gammaone/mu_bar;
					part.temp_ambient = temp0;		
					part.outflowFlag = (int) 2;
					j++;
				}
				else{
					for (k=0; k<ndim; k++) part.r[k] = rHotGas[ndim*j + k];
					for (k=0; k<ndim; k++) part.v[k] = 0.0;
					part.m = mp;
					part.u = temp0/gammaone/mu_bar;
					part.temp_ambient = temp0;		
					part.outflowFlag = (int) 0;
					j++;
				}
			}
			if(i >= Npart+NHotGas){
		  	for (k=0; k<ndim; k++) part.r[k] = r_shell[ndim*(i-(Npart+NHotGas)) + k];
		    for (k=0; k<ndim; k++) part.a[k] = (FLOAT) 0.0;	
				part.u = T_ambient/gammaone/mu_bar;
				part.temp_ambient = temp0;		
				part.m = mp;
				part.outflowFlag = (int) 2;
			}
		}
		
		sim->initial_h_provided = false;

		for (i=0; i<Npart+NHotGas-1; i++) {
			Particle<ndim>& part = hydro->GetParticlePointer(i);

        r_perp[0] = -part.r[1];
        r_perp[1] = part.r[0];
        r_perp[2] = 0.0;
        for (k=0; k<ndim; k++) part.v[k] = (r_perp[k]*omega*pc) / simunits.v.outscale;
      
		}

		delete[] v;
		delete[] r;
		delete[] rHotGas;
		delete[] r_shell;
		return;
	}

}


template class BonnorEbertSphereIc<1>;
template class BonnorEbertSphereIc<2>;
template class BonnorEbertSphereIc<3>;
