//=================================================================================================
//  AccretionFB.cpp
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



#include <iostream>
#include "AccretionFB.h"
#include "Constants.h"
#include "Hydrodynamics.h"
#include "Parameters.h"
#include "Particle.h"
#include "SimUnits.h"
#include "Nbody.h"
#include "Simulation.h"
using namespace std;

//=================================================================================================
//  AccretionFB::AccretionFB()
/// AccretionFB class constructor/destructor
//=================================================================================================
template <int ndim>
AccretionFB<ndim>::AccretionFB(SimUnits *simunits_, SimulationBase *sim_, string stellarEv_,
  const FLOAT gammam1_, const FLOAT mu_bar_, const FLOAT temp_ambient0_)
{
  DOUBLE num, denom, tempunit;
	FLOAT R_SUN;

  stellarEv = stellarEv_;
	simunits = simunits_;
	sim = sim_;


  // compute rad_const = Stefan-Boltzmann in code units
  num       = pow(simunits->r.outscale*simunits->r.outSI,2)*simunits->t.outscale*simunits->t.outSI;
  denom     = simunits->E.outscale * simunits->E.outSI;
  tempunit  = simunits->temp.outscale * simunits->temp.outSI;
  rad_const = stefboltz*(num*pow(tempunit,4.0))/denom;

  // compute G_Const in code units
  num       = pow(simunits->t.outscale*simunits->t.outSI,2)*simunits->m.outscale*simunits->m.outSI;
  denom     = pow(simunits->r.outscale*simunits->r.outSI,3);
  G_Const	  = G_const*num/denom;

  // compute L_sun in code units
  denom     = simunits->L.outscale * simunits->L.outSI;
  L_Sun     = L_sun/denom;

  // compute M_sun in code units
  denom     = simunits->m.outscale*simunits->m.outSI;
  M_Sun     = m_sun/denom;

  // compute R_sun in code units
  denom     = simunits->r.outscale*simunits->r.outSI;
  R_Sun     = r_sun/denom;

  // compute ambient temperature in code units
  temp_ambient0 	= temp_ambient0_/tempunit;

  gammam1   = gammam1_;
  mu_bar    = mu_bar_;

	cout << " R_Sun " << r_sun << " M_Sun " << m_sun << " L_Sun " << L_sun << " G_Const " << G_const << " rad_const " << stefboltz<< endl;
	cout << " R_Sun " << R_Sun << " M_Sun " << M_Sun << " L_Sun " << L_Sun << " G_Const " << G_Const << " rad_const " << rad_const<< endl;

}

template <int ndim>
AccretionFB<ndim>::~AccretionFB()
{
}

//=================================================================================================
//  AccretionFB::AmbientTemp()
/// Calculates ambient temperature for each particle
//=================================================================================================
template <int ndim>
void AccretionFB<ndim>::updateAccFBTemp(const int Npart, Hydrodynamics<ndim> *hydro, Sinks<ndim> *sinks)
{
    int Nsink = sinks->Nsink;                               			// Number of Sinks
		FLOAT T_max = 0.9*pow(10,6);
// Compute background temperature caused by stars and background radiation field.
//-------------------------------------------------------------------------------------------

#pragma omp parallel for default(none) shared(Nsink, sinks, hydro, cout, T_max)
    for (int i=0; i<Npart; i++){

    	int j;                                                          // Sink Particle cou$
    	int k;                                                          // Dimension counter
    	FLOAT T_bgr;                                                    // Summing up background tem$
    	FLOAT dist2;                                                    // Squared distance between $
    	FLOAT tempunit;                                                 // Unit scaling
    	FLOAT u_bgr;

    	Particle<ndim>& part = hydro->GetParticlePointer(i);


		for (j=0; j<Nsink; j++) {
			T_bgr = pow(temp_ambient0,4);
		    dist2 = 0.0;
		    for (k=0; k<ndim; k++) {
		            dist2 += pow(part.r[k] - sinks->sink[j].star->r[k],2);
		    }
				dist2 = max(dist2 , pow(R_Sun,2));
				T_bgr += sinks->sink[j].luminosity / (16* pi * rad_const * dist2);
		}
		T_bgr = pow(T_bgr,0.25);
		if (T_bgr > T_max) T_bgr = T_max;
		part.temp_ambient = T_bgr;

	}
    return;
}

//=================================================================================================
//  EpisodicAFB::ContinuousAFB()
/// EpisodicAFB class constructor/destructor
//=================================================================================================
template <int ndim>
ContinuousAFB<ndim>::ContinuousAFB(SimUnits *simunits_, SimulationBase *sim_, string stellarEv_, const FLOAT gammam1_, const FLOAT mu_bar_, const FLOAT temp_ambient0_,
	const FLOAT dmdt_max_, const FLOAT f_angmom_, const FLOAT AccFBminMass_):
	AccretionFB<ndim>(simunits_, sim_, stellarEv_, gammam1_, mu_bar_, temp_ambient0_)
{
	AccFBminMass = AccFBminMass_;
	dmdt_max = dmdt_max_;
	f_angmom = f_angmom_;
}

template <int ndim>
ContinuousAFB<ndim>::~ContinuousAFB()
{
}

//=================================================================================================
//  AccretionFB::Luminosity()
/// Calculates the Luminosity of each star
//=================================================================================================
template <int ndim>
void ContinuousAFB<ndim>::Luminosity(Hydrodynamics<ndim> *hydro, Sinks<ndim> *sinks, const FLOAT f_acc)
{
	int Nsink = sinks->Nsink;			// Number of sink particle
	int i;								// Sink particle counter

// Compute accretion rate
//-------------------------------------------------------------------------------------------

	Accretion(sinks, hydro);

// Compute Luminosity
//-------------------------------------------------------------------------------------------
  if (stellarEv == "none"){

#pragma omp parallel for default(none) shared(Nsink, sinks) private(i)
  	for (i=0; i<Nsink; i++) {
  		if (sinks->sink[i].star->m < AccFBminMass){
  			sinks->sink[i].acc_flag = 0.0;
  			continue;
  		}
  		sinks->sink[i].luminosity = pow(sinks->sink[i].M_star/M_Sun,3)*L_Sun+f_acc*(G_Const*sinks->sink[i].M_star_old*
  			sinks->sink[i].dmdt_star)/(3*R_Sun);
  	}
  }
return;
}

//=================================================================================================
//  AccretionFB::Accretion()
/// Smoothes the accretion rate
//=================================================================================================
template <int ndim>
void ContinuousAFB<ndim>::Accretion(Sinks<ndim> *sinks, Hydrodynamics<ndim> *hydro)
{
	int i;													// Sink particle counter
	int Nsink = sinks->Nsink;				// Number of Sinks
	FLOAT t = sim->t;
	FLOAT dm;

	// compute mass of a single SPH-particle in code units
	Particle<ndim>& part = hydro->GetParticlePointer(1);
	dm = part.m;

#pragma omp parallel for default(none) shared(Nsink, sinks, t, dm) private(i)
	for (i=0; i<Nsink; i++) {



		FLOAT dt;
		FLOAT M_IAD;
		FLOAT dmdt;
		FLOAT dmdt_IAD;

		dt = t - sinks->sink[i].t_old;
		M_IAD = sinks->sink[i].star->m - sinks->sink[i].M_star;
		dmdt = M_IAD/dt;

		sinks->sink[i].m_old = sinks->sink[i].star->m;
		sinks->sink[i].M_star_old = sinks->sink[i].M_star;
		sinks->sink[i].M_star = sinks->sink[i].star->m;
		sinks->sink[i].dmdt_star = dmdt;

  	// Angular Momentum to be carried in Outflows (todo: turn off if Outflows are turned of!)
		//sinks->sink[i].dL_Outflow = f_angmom * sqrt(pow(sinks->sink[i].angmom[0] - sinks->sink[i].L_star[0],2)
		//	+ pow(sinks->sink[i].angmom[1] - sinks->sink[i].L_star[1],2)
		//	+ pow(sinks->sink[i].angmom[2] - sinks->sink[i].L_star[2],2))/ max(1.0,floor(M_IAD/dm));
		//Update sink properties for next step
		sinks->sink[i].acc_flag = 1;
		sinks->sink[i].t_old = t;
	}
		//cout << sinks->sink[0].dL_Outflow << endl;
	return;
}

//=================================================================================================
//  EpisodicAFB::EpisodicAFB()
/// EpisodicAFB class constructor/destructor
//=================================================================================================
template <int ndim>
EpisodicAFB<ndim>::EpisodicAFB(SimUnits *simunits_, SimulationBase *sim_, string stellarEv_, const FLOAT AccFBminMass_,
	const FLOAT dmdt_BRG_, const FLOAT alpha_MRI_, const FLOAT gammam1_, const FLOAT mu_bar_, const FLOAT temp_ambient0_, const FLOAT f_angmom_):
	AccretionFB<ndim>(simunits_, sim_, stellarEv_, gammam1_, mu_bar_, temp_ambient0_)
{
  FLOAT num, denom;

	AccFBminMass = AccFBminMass_;
	alpha_MRI = alpha_MRI_;
	f_angmom = f_angmom_;

// Compute required scaling quantities
//-------------------------------------------------------------------------------------------


  // compute quiescent accretion rate in code units
  denom     = simunits->dmdt.outscale;
  dmdt_BRG 	= (dmdt_BRG_)/denom;

  // compute time unit
	yr = 1.0e-6/simunits->t.outscale;
	dmdt_MRI = (alpha_MRI/0.1)*(M_Sun/yr)*5.0e-4;

}

template <int ndim>
EpisodicAFB<ndim>::~EpisodicAFB()
{
}

//=================================================================================================
//  AccretionFB::Luminosity()
/// Calculates the Luminosity of each star
//=================================================================================================
template <int ndim>
void EpisodicAFB<ndim>::Luminosity(Hydrodynamics<ndim> *hydro, Sinks<ndim> *sinks_, const FLOAT f_acc)
{
	sinks = sinks_;
	int Nsink = sinks->Nsink;			// Number of sink particle
	int i;								// Sink particle counter

// Call Eppisodic Accretion module
//-------------------------------------------------------------------------------------------
	Accretion(sinks, hydro);

// Compute Luminosity
//-------------------------------------------------------------------------------------------

  if (stellarEv != "none"){
    return;
  }
#pragma omp parallel for default(none) shared(Nsink) private(i)
  for (i=0; i<Nsink; i++) {
    if (sinks->sink[i].star->m < AccFBminMass){
      sinks->sink[i].acc_flag = 0.0;
      continue;
    }
  sinks->sink[i].luminosity = pow(sinks->sink[i].M_star/M_Sun,3)*L_Sun +
    f_acc*(G_Const*sinks->sink[i].M_star*sinks->sink[i].dmdt_star)/(3*R_Sun);
  }

  return;
}



//=================================================================================================
//  EpisodicAFB::Accretion()
/// Calculates episodic accretion rate
//=================================================================================================
template <int ndim>
void EpisodicAFB<ndim>::Accretion(Sinks<ndim> *sinks, Hydrodynamics<ndim> *hydro)
{
	int i;																	// Sink Particle counter
	int Nsink = sinks->Nsink;								// Number of Sinks
	FLOAT t	 	= sim->t;											// Simulation time
	FLOAT dm;																// mass of a single SPH particle

	// compute mass of a single SPH-particle in code units
	Particle<ndim>& part = hydro->GetParticlePointer(1);
	dm = part.m;

#pragma omp parallel for default(none) shared(Nsink, sinks, hydro, cout, t, dm) private(i)
	for (i=0; i<Nsink; i++) {

		FLOAT M_MRI; 													// Mass of the outburst
		FLOAT t_MRI;													// timescale of an outburst
		FLOAT t0;															// Simulation time of the start of the last outburst

		FLOAT dmdt;														// Accretionrate onto the sink
		FLOAT dt;															// Timestep of the sink
		FLOAT dmdt_max = 1.0e8*(M_Sun/yr);		// Max accretionrate
		FLOAT dmdt_EA;												// Average accretionrate during an outburst
		FLOAT dmdt_star;											// Accretionrate onto the inner star
		FLOAT M_IAD;													// Mass of the inner accretion disc
    FLOAT alpha;                          // fraction of accreted mass

		//Update sink variables
		t0 = sinks->sink[i].t0;
    dt =  sinks->sink[i].dt;
		dmdt = (sinks->sink[i].star->m - sinks->sink[i].m_old)/dt; // dt_old
		dmdt_EA = 0.0;
		sinks->sink[i].M_star_old = sinks->sink[i].M_star;

		if (sinks->sink[i].star->m < AccFBminMass){
			sinks->sink[i].M_star = sinks->sink[i].star->m;
      for (int k=0; k<ndim; k++) sinks->sink[i].L_star[k] = sinks->sink[i].angmom[k];
			sinks->sink[i].acc_flag = 0.0;
			continue;
		}

		//Calculate Outburst properties: Mass of the inner accretiondisc (M_IAD) and timescale of the outburst (t_MRI)
		M_IAD = sinks->sink[i].m_old - sinks->sink[i].M_star;
		t_MRI = (0.25*1000)*yr*(0.1/alpha_MRI)*pow((sinks->sink[i].M_star/0.2*M_Sun),2.0/3.0)
			*pow(dmdt*yr/(M_Sun*10e-5),1.0/9.0);
		M_MRI = dmdt_MRI * t_MRI;


//-------------------------------------------------------------------------------------------
// If there is no outburst atm.:
//-------------------------------------------------------------------------------------------

		if (sinks->sink[i].acc_flag == 0){

// Check for outburst
//-------------------------------------------------------------------------------------------

			if (M_IAD > M_MRI and sinks->sink[i].star->m > AccFBminMass and dmdt > 0 and M_MRI > 0){
				sinks->sink[i].acc_flag = 1;
				sinks->sink[i].t0 = t;
				sinks->sink[i].DT_MRI = t_MRI;

  			// Angular Momentum to be carried in Outflows (todo: turn off if Outflows are turned of!)
				//sinks->sink[i].dL_Outflow = f_angmom * sqrt(pow(sinks->sink[i].angmom[0] - sinks->sink[i].L_star[0],2)
				//	+ pow(sinks->sink[i].angmom[1] - sinks->sink[i].L_star[1],2)
				//	+ pow(sinks->sink[i].angmom[2] - sinks->sink[i].L_star[2],2))/ floor(M_MRI/dm);


				if (1.58*dmdt_MRI < dmdt_max){
					sinks->sink[i].dmdt_MRI_0 = 1.58*dmdt_MRI;
				}
				else{
					sinks->sink[i].dmdt_MRI_0 = dmdt_max;
				}

				dmdt_EA = sinks->sink[i].dmdt_MRI_0;


			}
			dmdt_star = dmdt_BRG + dmdt_EA;
			sinks->sink[i].dmdt_star = dmdt_star;

// update M_IAD
//-------------------------------------------------------------------------------------------
      alpha = (sinks->sink[i].dmdt_star-dmdt)*dt / M_IAD;
			M_IAD -= alpha*M_IAD;

			if(M_IAD < 0.0){
				M_IAD = 0.0;								// Check if acc_rate to high if MIAD < 0
				cout << "MIAD < 0.0 M_sun " << endl;
			}
			if(t-t0 >= sinks->sink[i].DT_MRI or M_IAD <= 0.0){
				sinks->sink[i].acc_flag = 0;
			}
      else{
		    sinks->sink[i].M_star = sinks->sink[i].star->m - M_IAD;
        for (int k=0; k<ndim; k++) sinks->sink[i].L_star[k] += alpha*sinks->sink[i].L_IAD[k];
        for (int k=0; k<ndim; k++) sinks->sink[i].L_IAD[k] -= alpha*sinks->sink[i].L_IAD[k];
      }
  }



//-------------------------------------------------------------------------------------------
// If there is an outburst atm.:
//-------------------------------------------------------------------------------------------
		else if (sinks->sink[i].acc_flag == 1){

// check for restart of the outburst
//-------------------------------------------------------------------------------------------

			if (M_IAD > M_MRI and sinks->sink[i].star->m > AccFBminMass and dmdt > 0 and M_MRI > 0){

  			// Angular Momentum to be carried in Outflows (todo: turn off if Outflows are turned of!)
				//sinks->sink[i].dL_Outflow = f_angmom * sqrt(pow(sinks->sink[i].angmom[0] - sinks->sink[i].L_star[0],2)
				//	+ pow(sinks->sink[i].angmom[1] - sinks->sink[i].L_star[1],2)
				//	+ pow(sinks->sink[i].angmom[2] - sinks->sink[i].L_star[2],2))/ floor(M_MRI/dm);

				sinks->sink[i].acc_flag = 1;
				sinks->sink[i].t0 = t;
				sinks->sink[i].DT_MRI = t_MRI;

				if (1.58*dmdt_MRI < dmdt_max){
					sinks->sink[i].dmdt_MRI_0 = 1.58*dmdt_MRI;
				}
				else{
					sinks->sink[i].dmdt_MRI_0 = dmdt_max;
				}
			}
			if ((t-t0) < sinks->sink[i].DT_MRI) {
				dmdt_EA = sinks->sink[i].dmdt_MRI_0 *exp(-(t-t0)/sinks->sink[i].DT_MRI);
			}
			else {
				dmdt_EA = 0.0;
			}
			dmdt_star = dmdt_BRG + dmdt_EA;
			sinks->sink[i].dmdt_star = dmdt_star;

      alpha = (sinks->sink[i].dmdt_star-dmdt)*dt / M_IAD;
			M_IAD -= alpha*M_IAD;

			if(M_IAD < 0.0){
				M_IAD = 0.0;								// Check if acc_rate to high if MIAD < 0
				cout << "MIAD < 0.0 M_sun " << endl;
			}
			if(t-t0 >= sinks->sink[i].DT_MRI or M_IAD <= 0.0){
				sinks->sink[i].acc_flag = 0;
			}
      else{
		    sinks->sink[i].M_star = sinks->sink[i].star->m - M_IAD;
        for (int k=0; k<ndim; k++) sinks->sink[i].L_star[k] += alpha*sinks->sink[i].L_IAD[k];
        for (int k=0; k<ndim; k++) sinks->sink[i].L_IAD[k] -= alpha*sinks->sink[i].L_IAD[k];
      }
		}

	//Update sink properties for next step
	sinks->sink[i].t_old = t;
	sinks->sink[i].m_old = sinks->sink[i].star->m;
	}

	return;
}

template class AccretionFB<1>;
template class AccretionFB<2>;
template class AccretionFB<3>;
template class ContinuousAFB<1>;
template class ContinuousAFB<2>;
template class ContinuousAFB<3>;
template class EpisodicAFB<1>;
template class EpisodicAFB<2>;
template class EpisodicAFB<3>;


