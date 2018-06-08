//=================================================================================================
//  EnergyRadws.cpp
//  Contains functions for Rad-WS radiation cooling scheme (Stamatellos et al. 2007) class.
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


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

#include "../Headers/Integration.h"
#include "Constants.h"
#include "Debug.h"
#include "EnergyEquation.h"
#include "EOS.h"
#include "Exception.h"
#include "Hydrodynamics.h"
#include "Parameters.h"
#include "Particle.h"
#include "SimUnits.h"
#include "SmoothingKernel.h"
using namespace std;

FLOAT bdens;
FLOAT btemp;
FLOAT densmin;
FLOAT densmax;
FLOAT tempmin;
FLOAT tempmax;

FLOAT r_unit      ;
FLOAT time_unit   ;
FLOAT dens_unit   ;
FLOAT temp_unit   ;
FLOAT energy_unit ;
FLOAT u_unit      ;
FLOAT kappa_unit  ;

//=================================================================================================
//  EnergyRadws::EnergyRadws()
/// EnergyRadws class constructor
//=================================================================================================
template <int ndim, template <int> class ParticleType>
EnergyRadws<ndim,ParticleType>::EnergyRadws
 ( DOUBLE energy_mult_, string radws_table, FLOAT temp_ambient_,
  SimUnits *simunits, EOS<ndim> *eos_) :
  EnergyEquation<ndim>(energy_mult_)
{
  int i, j, l;
  // int ndens, ntemp; defined in EnergyEquation
  DOUBLE eos_dens_, eos_temp_, eos_energy_, eos_mu_, kappa_, kappar_, kappap_, eos_gamma_;
  DOUBLE num, denom, tempunit;
  string line;

  eos = eos_;

  // compute rad_const = Stefan-Boltzmann in code units
  num       = pow(simunits->r.outscale*simunits->r.outSI,2)*simunits->t.outscale*simunits->t.outSI;
  denom     = simunits->E.outscale * simunits->E.outSI;
  tempunit  = simunits->temp.outscale * simunits->temp.outSI;
  rad_const = stefboltz*(num*pow(tempunit,4.0))/denom;

  dens_unit = simunits->rho.outscale * simunits-> rho.outcgs;
  temp_unit = simunits->temp.outscale * simunits->temp.outcgs;
  u_unit = simunits->u.outscale * simunits-> u.outcgs;

  temp_ambient0 = temp_ambient_ / tempunit;
  cout << "\n\n TEMP UNIT = " << tempunit << "\n\n";
  cout << "Temp_ambient : " << temp_ambient0 << endl;
  cout << "units     : " << num << " "<< denom << " " << tempunit << endl;
  cout << "u_unit    : " << simunits->u.outscale * simunits->u.outcgs << endl;
  cout << "time_unit : " << simunits->t.outscale * simunits->t.outcgs << endl;
  cout << "rad_const : " << rad_const << endl;
  cout << "stefboltz : " << stefboltz << endl;
  cout << "rad_constscaled : " << rad_const*(simunits->E.outscale*simunits->E.outSI/
    (pow(simunits->r.outscale*simunits->r.outSI,2)*simunits->t.outscale*simunits->t.outSI*
     pow(simunits->temp.outscale*simunits->temp.outSI,4))) << endl;

  // check that user wants cm2_g for opacity units
  if (simunits ->kappa.outunit != "cm2_g") {
    cout << "ERROR! Selected wrong unit for opacity. Use cm2_g" <<endl;
    cout << simunits->kappa.outunit << endl;
    ExceptionHandler::getIstance().raise("Error : Wrong units for opacity.  Use cm2_g");
  }

  ifstream file;
  file.open(radws_table.c_str(), ios::in);
  // allocate table entries


  cout << "EnergyRadws : " << radws_table << "  "<< file.good() << endl;

  //-----------------------------------------------------------------------------------------------
  if (file.good()) {
    getline(file, line);
    cout << line << endl;
    istringstream istr(line);
    istr >> ndens >> ntemp ;

    cout << "ndens, ntemp = " << ndens << "  " << ntemp << endl;

    eos_dens     = new FLOAT[ndens];
    eos_temp     = new FLOAT[ntemp];
    eos_energy   = new FLOAT*[ndens];
    eos_mu       = new FLOAT*[ndens];
    eos_gamma    = new FLOAT*[ndens];
    kappa_table  = new FLOAT*[ndens];
    kappar_table = new FLOAT*[ndens];
    kappap_table = new FLOAT*[ndens];

    for (i=0; i<ndens; i++){
      eos_energy[i]   = new FLOAT[ntemp];
      eos_mu[i]       = new FLOAT[ntemp];
      eos_gamma[i]    = new FLOAT[ntemp];
      kappa_table[i]  = new FLOAT[ntemp];
      kappar_table[i] = new FLOAT[ntemp];
      kappap_table[i] = new FLOAT[ntemp];
    }

    // read table
    i = 0;
    l = 0;
    j = 0;

    //---------------------------------------------------------------------------------------------
    while (getline(file, line)) {
      istringstream istr(line);

      if (istr >> eos_dens_ >> eos_temp_ >> eos_energy_ >> eos_mu_>> kappa_ >> kappar_ >> kappap_ >> eos_gamma_) {

        eos_energy[i][j]   = eos_energy_/(simunits->u.outscale * simunits-> u.outcgs);
        eos_mu[i][j]       = eos_mu_;
        eos_gamma[i][j]    = eos_gamma_;
        kappa_table[i][j]  = kappa_/(simunits->kappa.outscale * simunits-> kappa.outcgs);
        kappar_table[i][j] = kappar_/(simunits->kappa.outscale * simunits-> kappa.outcgs);
        kappap_table[i][j] = kappap_/(simunits->kappa.outscale * simunits-> kappa.outcgs);

        if (l < ntemp) {
          eos_temp[l] = log10(eos_temp_/(simunits->temp.outscale * simunits->temp.outcgs));
        }

        ++l;
        ++j;

        if (!(l % ntemp)) {
          eos_dens[i] = log10(eos_dens_/(simunits->rho.outscale * simunits-> rho.outcgs));
          ++i;
          j = 0;
        }

      } else cout << "Dateifehler" << endl;

    }
    //---------------------------------------------------------------------------------------------

    file.close();
  }
  //-----------------------------------------------------------------------------------------------

  cout << "eos_dens = " << eos_dens[0] << " - \t" << eos_dens[ndens - 1] << endl;
  cout << "eos_temp = " << eos_temp[0] << " - \t" << eos_temp[ntemp - 1] << endl;
}



//=================================================================================================
//  EnergyRadws::~EnergyRadws()
/// EnergyRadws class destructor
//=================================================================================================
template <int ndim, template <int> class ParticleType>
EnergyRadws<ndim,ParticleType>::~EnergyRadws()
{
  for (int i = 0; i < ndens; i++) {
    delete[] eos_energy[i];
    delete[] eos_mu[i];
    delete[] kappa_table[i];
    delete[] kappar_table[i];
    delete[] kappap_table[i];
  }

  delete[] eos_energy;
  delete[] eos_mu;
  delete[] kappa_table;
  delete[] kappar_table;
  delete[] kappap_table;

}



//=================================================================================================
//  EnergyRadws::EnergyIntegration
/// Integrate internal energy to first order from the beginning of the step to
/// the current simulation time, i.e. u(t+dt) = u(t) + dudt(t)*dt .
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyRadws<ndim,ParticleType>::EnergyIntegration
 (const int n,                         ///< [in] Integer time in block time struct
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  Hydrodynamics<ndim>* hydro)
{
  //cout << "int1" << endl;
  int i;                               // Particle counter
  FLOAT dt;                            // Timestep since start of step
  //ParticleType<ndim>* partdata = hydro->template GetParticleArray<ParticleType>();
  //MeshlessFV<ndim>* mfv = reinterpret_cast<MeshlessFV<ndim>*>(hydro);
  //MeshlessFVParticle<ndim>* partdata = mfv->GetMeshlessFVParticleArray()
  //debug2("[EnergyRadws::EnergyIntegration]");
  //CodeTiming::BlockTimer timer = timing->StartNewTimer("ENERGY_RADWS_INTEGRATION");

  //cout << "int2" << endl;

  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dt,i) shared(hydro,cout) //partdata,
  for (i=0; i<hydro->Nhydro; i++) {
    //ParticleType<ndim>& part = partdata[i];
  	//cout << "int3" << endl;
  	//FLOAT u0;
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    if (part.flags.is_dead()) continue;

  	//cout << "int4" << endl;
    // Compute time since beginning of current step
    //dt = t - part.tlast;
    int dn = n - part.nlast;
    FLOAT dt = timestep*(FLOAT) dn;

    //cout << dt << " " << part.dt << " " << part.dt_next << " " << n << " " << part.nstep << " " << part.nlast <<  endl;
    //dt = part.dt;

	//if (dn <= part.nstep){
    //u0 = part.u;
    if (part.dt_therm <= small_number) {
      part.u = part.u0; // MERCER
    }
    else if (dt < (FLOAT) 40.0 * part.dt_therm) {
      part.u = part.u * exp(-dt / part.dt_therm)
             + part.ueq * ((FLOAT) 1.0 - exp(-dt / part.dt_therm));
    }
    else if (dt >= (FLOAT) 40.0 * part.dt_therm) {
      part.u = part.ueq;
    }
    GetTemp(part);
    //cout  << " u " << part.u << " ueq " << part.ueq << " dt_therm " << part.dt_therm << " dudt " << part.dudt <<  endl;
    //part.du_rad = part.u - u0;
    //cout << " part.du_rad " << part.du_rad << endl;
    //if (part.ueq >100000) cout << part.ueq << endl;
	//}
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  EnergyRadws::EndTimestep
/// Record all important thermal quantities at the end of the step for start of the new timestep.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyRadws<ndim,ParticleType>::EndTimestep
 (const int n,                         ///< [in] Integer time in block time struct
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  Hydrodynamics<ndim>* hydro)
{
  int dn;                              // Integer time since beginning of step
  int i;                               // Particle counter
  FLOAT temp;                          // ..
  FLOAT temp_ambient;

  //FLOAT dt_therm;                      // ..
  //ParticleType<ndim>* partdata = hydro->template GetParticleArray<ParticleType>();
  //cout << "end1" << endl;
  //debug2("[EnergyRadws::EndTimestep]");
  //CodeTiming::BlockTimer timer = timing->StartNewTimer("ENERGY_RADWS_END_TIMESTEP");

  //cout << "end" << endl;
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,temp,temp_ambient) shared(hydro,cout,temp_unit) // partdata,
  for (i=0; i<hydro->Nhydro; i++) {
    //ParticleType<ndim> &part = partdata[i];
  //cout << "end2" << endl;
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    if (part.flags.is_dead()) continue;
    dn = n - part.nlast;
  ///cout << "end3" << endl;
    FLOAT r_centre;
    FLOAT gpot;

    if (part.flags.check(end_timestep)) {

    	temp_ambient =  max(part.temp_ambient,temp_ambient0);
		gpot = part.gpot;
		/*if (sinks->Nsink > 0){

	  	for (int j=0; j<sinks->Nsink; j++) {
			r_centre = 0.0;
	  		for (int k=0; k<3; k++) r_centre += pow(sinks->sink[j].star->r[k]-part.r[k],2);
				r_centre = sqrt(r_centre);
				if (r_centre < sinks->sink[j].radius) r_centre = sinks->sink[j].radius;
				//cout << " gpot 1 : " << gpot << endl;
				gpot -= sinks->sink[j].star->m/r_centre; // minus / plus ?
				if (gpot < 1.0)	gpot = 1.0;
				//cout << " gpot 2 : " << gpot << " " << sinks->sink[j].star->m << " " << r_centre <<  endl;
			}
		}*/


      // MERCER
      // Instead of doing this, find temperature in the table from part.u and part.rho.
      // temp = eos->Temperature(part);
      temp = GetTemp(part); // MERCER


      // Get new ueq and dt_therm
     //cout << temp << " " << part.ueq << " " << part.u << " " << part.dudt <<  endl;
      EnergyFindEqui(part.rho, temp, gpot, part.u, temp_ambient, part.dudt,
                     part.ueq, part.dt_therm);


      part.u0 = part.u;
      part.dudt0 = part.dudt; // MERCER
    }
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  EnergyRadws::~EnergyFindEqui()
/// Computes the thermal equilibrium state of the particle (i.e. its equilibrium temperature and
/// internal energy) including the thermal timescale to reach this equilibrium.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyRadws<ndim,ParticleType>:: EnergyFindEqui
 (const FLOAT rho,                             ///< [in] ..
  const FLOAT temp,                            ///< [in] ..
  const FLOAT gpot,                            ///< [in] ..
  const FLOAT u,                               ///< [in] ..
	const FLOAT temp_ambient,										 ///< [in] ..
  FLOAT dudt,                            			 ///< [in] ..
  FLOAT &ueq_p,                                ///< [out] ..
  FLOAT &dt_thermal)                           ///< [out] ..
{
  const FLOAT fcolumn2 = 0.010816;             // ..
  const FLOAT col2     = fcolumn2*gpot*rho;    // ..
  const FLOAT logrho   = log10(rho);           // ..
  const int idens      = GetIDens(logrho);     // ..
  int itemp;                                   // ..
  FLOAT logtemp;                               // ..
  FLOAT Tequi;                                 // ..
  FLOAT dudt_eq;                               // ..
  FLOAT dudt_rad;                              // ..
  FLOAT dudt_tot;                              // ..
  FLOAT kappa;                                 // ..
  FLOAT kappar;                                // ..
  FLOAT kappap;                                // ..

  logtemp = log10(temp);
  itemp   = GetITemp(logtemp);

  assert(idens >= 0 && idens <= ndens - 1);
  assert(itemp >= 0 && itemp <= ntemp - 1);

  GetKappa(idens, itemp, logrho, logtemp, kappa, kappar, kappap);
  dudt_rad = ebalance((FLOAT) 0.0, temp_ambient, temp, kappa, kappap, col2);

  // Calculate equilibrium temperature using implicit scheme described in Stamatellos et al. (2007)
  EnergyFindEquiTemp(idens, rho, temp, col2, dudt, temp_ambient, Tequi);
  assert(Tequi >= temp_ambient);

  // Get ueq_p and dudt_eq from Tequi
  logtemp = log10(Tequi);
  itemp = GetITemp(logtemp);
  GetKappa(idens, itemp, logrho, logtemp, kappa, kappar, kappap);

  ueq_p = GetEnergy(idens, itemp, logrho, logtemp); // MERCER
  dudt_eq = ebalance((FLOAT) 0.0 ,temp_ambient ,Tequi, kappa, kappap, col2);

  // Thermalization time scale
  dudt_tot = -(dudt_eq - dudt_rad);
  if (dudt_tot == 0.0) {
    dt_thermal = 1.e30;
  }
  else {
    dt_thermal = (ueq_p - u) / (dudt + dudt_rad);
  }

  //cout << " dt_thermal " << dt_thermal << " temp " << temp << " u " << u << " ueq_p " << ueq_p << " dudt " << dudt << " temp_ambient " << temp_ambient0 << " dudt_rad " << dudt_rad << " idens " << idens <<  " itemp " << itemp << endl;
  return;
}



//=================================================================================================
//  EnergyRadws::EnergyFindEquiTemp
/// EnergyFindEquiTemp returns equilibrium temperature
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyRadws<ndim,ParticleType>::EnergyFindEquiTemp
 (const int idens,                         ///< ..
  const FLOAT rho,                         ///< ..
  const FLOAT temp,                        ///< ..
  const FLOAT col2,                        ///< ..
  const FLOAT dudt,                        ///< ..
  const FLOAT temp_ambient,                ///< ..
  FLOAT &Tequi)                            ///< ..
                      ///< ..
{
  int itemp;
  int itemplow, itemphigh;

  FLOAT accuracy = 0.001; // not shure what is a good accuracy (Paul)
  FLOAT balance, minbalance;
  FLOAT balanceLow, balanceHigh ;
  FLOAT kappa, kappar, kappap;
  FLOAT kappaLow, kappaHigh, kappapLow, kappapHigh;
  FLOAT logtemp, logrho;
  FLOAT Tlow, Thigh, Tlow_log, Thigh_log, Tequi_log;
  FLOAT mu_bar_high, mu_bar_low;
  FLOAT dtemp;

  logrho = log10(rho);

  // Get min. kappa and min balance
  logtemp = log10(temp_ambient);
  itemp = GetITemp(logtemp);
  GetKappa(idens, itemp, logrho, logtemp, kappa, kappar, kappap);
  minbalance = ebalance((FLOAT) 0.0, temp_ambient, temp_ambient, kappa, kappap, col2);

  logtemp = log10(temp);
  itemp = GetITemp(logtemp);
  GetKappa(idens, itemp, logrho, logtemp, kappa, kappar, kappap);
  balance = ebalance((FLOAT) 0.0, temp_ambient, temp, kappa, kappap, col2);



  // Find equilibrium temperature
  //cout << dudt << " " << -balance << endl;
  // ------------------------------------------------------------------------------------------------
  if (dudt < -balance) {
  	//cout <<  dudt << " " <<   -balance << endl;            // <=
    if (dudt <= -minbalance){
    //if (dudt <= 0.0){
      Tequi = temp_ambient;
      return;
    }

    // COOLING
    if (itemp <= 1) {
      cout << "Leaving immediately! // temp = " << temp << endl;
      Tequi = pow(10.0, eos_temp[1]);
      return;
    }

    Tlow_log  = eos_temp[itemp];
    Tlow      = pow(10.0, Tlow_log);
    Thigh_log = eos_temp[itemp + 1];
    Thigh     = pow(10.0, Thigh_log);

    itemplow  = itemp;
    GetKappa(idens, itemplow, logrho, Tlow_log, kappaLow, kappar, kappapLow);
    balanceLow = ebalance(dudt, temp_ambient, Tlow, kappaLow, kappapLow, col2);

    itemphigh = itemp+1;
    GetKappa(idens, itemphigh, logrho, Thigh_log, kappaHigh, kappar, kappapHigh);
    balanceHigh = ebalance(dudt, temp_ambient ,Thigh, kappaHigh, kappapHigh, col2);


    while (balanceLow * balanceHigh > 0.0){

      if (itemplow <= 1) { // MERCER itemp to itemplow
        cout << "Reached itemp=1, returning " << dudt << endl;
        Tequi = pow(10.0, temp_ambient);
        return;
      }

      Thigh       = Tlow;
      kappaHigh   = kappaLow;
      kappapHigh  = kappapLow;
      balanceHigh = balanceLow;
      itemplow    = itemplow - 1;
      Tlow_log    = eos_temp[itemplow];
      Tlow        = pow(10.0, Tlow_log);
		//cout << balanceLow <<  " " << balanceHigh << " " << temp << " " << itemp << " " << itemplow << " " << itemphigh << " " << dudt << " " << Tlow  <<  endl;
      GetKappa(idens, itemplow, logrho, Tlow_log, kappaLow, kappar, kappapLow);
      balanceLow = ebalance(dudt, temp_ambient , Tlow, kappaLow, kappapLow, col2);
    }

    itemphigh = itemplow + 1;
  } else {
    // HEATING: ERROR MAY BE HIDING HERE
    if (itemp >= ntemp - 2) {
      Tequi = pow(10.0, eos_temp[ntemp - 2]);
      return;
    }

    Tlow_log  = eos_temp[itemp -1];
    Tlow      = pow(10.0, Tlow_log);
    Thigh_log = eos_temp[itemp];
    Thigh     = pow(10.0, Thigh_log);

    itemplow = itemp - 1;
    GetKappa(idens, itemplow, logrho, Tlow_log, kappaLow, kappar, kappapLow);
    balanceLow = ebalance(dudt, temp_ambient, Tlow, kappaLow, kappapLow, col2);

    itemphigh = itemp;
    GetKappa(idens, itemphigh, logrho, Thigh_log, kappaHigh, kappar, kappapHigh);
    balanceHigh = ebalance(dudt, temp_ambient, Thigh, kappaHigh, kappapHigh, col2);

    while (balanceLow * balanceHigh > 0.0) {
      if (itemphigh >= ntemp - 2) { // MERCER itemp to itemphigh
        cout << "itemp reached max " << itemphigh << " " << itemp << " " << temp << endl;
        cout << balanceLow <<  " " << balanceHigh << " " << balance << " " <<  temp << " " << itemp << " " << itemplow << " " << itemphigh << " " << dudt << " " << Tlow << " " << Thigh << endl;
        Tequi = pow(10.0, eos_temp[ntemp - 2]);
        return;
      }

      Tlow       = Thigh;
      kappaLow   = kappaHigh;
      kappapLow  = kappapHigh;
      balanceLow = balanceHigh;
      itemphigh  = itemphigh + 1;
      Thigh_log  = eos_temp[itemphigh];
      Thigh      = pow(10.0, Thigh_log);

      GetKappa(idens, itemphigh, logrho, Thigh_log, kappaHigh, kappar, kappapHigh);
      balanceHigh = ebalance(dudt, temp_ambient, Thigh, kappaHigh, kappapHigh, col2);
    }

    itemplow = itemphigh - 1;
  }

  Tequi = 0.5 * (Tlow + Thigh);
  dtemp = Thigh - Tlow;

  // Refine the search in between Thigh and Tlow
  //-----------------------------------------------------------------------------------------------
  while (dtemp != 0.0 && fabs(2.0 * dtemp / (Thigh + Tlow)) > accuracy) {
    Tequi_log = log10(Tequi);
    GetKappa(idens, itemplow, logrho, Tlow_log, kappaLow, kappar, kappapLow);
    balance = ebalance(dudt, temp_ambient , Tequi, kappa, kappap, col2);

    if (balance == 0.0) {
      Tequi = 0.5 * (Thigh + Tlow);
      return;
    }

    if (balanceLow * balance < 0.0){
      Thigh       = Tequi;
      balanceHigh = balance;
      kappaHigh   = kappa;
      kappapHigh  = kappap;
    }
    else {
      Tlow       = Tequi;
      balanceLow = balance;
      kappaLow   = kappa;
      kappapLow  = kappap;
    }

    Tequi = 0.5 * (Thigh + Tlow);
    //mu_bar_equi = (mu_bar_high + mu_bar_low) / 2;
    kappa = (kappaHigh + kappaLow) / 2;
    kappap = (kappapHigh + kappapLow) / 2;
    dtemp = Thigh - Tlow;
  }

  if (Tequi < temp_ambient) {
    Tequi = temp_ambient;
  }

  Tequi_log   = log10(Tequi);
  //mu_bar_equi = GetMuBar(idens, itemplow, logrho, Tequi_log);

  return;
}



// ------------------------------------------------------------------------------------------//
//			find  closest index in list for value: level
// ------------------------------------------------------------------------------------------//

template <typename BidirectionalIterator, typename T>
BidirectionalIterator getClosest(BidirectionalIterator first,
                                 BidirectionalIterator last,
                                 const T &value)
{
  BidirectionalIterator before = std::lower_bound(first, last, value);

  if (before == first) return first;
  if (before == last)  return --last; // iterator must be bidirectional

  BidirectionalIterator after = before;
  --before;

  return (*after - value) < (value - *before) ? after : before;
}


template <typename BidirectionalIterator, typename T>
std::size_t getClosestIndex(BidirectionalIterator first,
                            BidirectionalIterator last,
                            const T &value)
{
  return std::distance(first, getClosest(first, last, value));
}



//=================================================================================================
//  EnergyRadws::GetIDens()
/// GetIDens returns table index for density
//=================================================================================================
template <int ndim, template <int> class ParticleType>
int EnergyRadws<ndim,ParticleType>::GetIDens
 (const FLOAT rho)
{
  int idens = getClosestIndex(eos_dens, eos_dens + ndens, rho);

  if (rho < eos_dens[idens]) {
    idens = max(idens - 1, 0);
  }
  return idens;
}



//=================================================================================================
//  EnergyRadws::GetITemp()
/// GetITemp returns table index for temperature
//=================================================================================================
template <int ndim, template <int> class ParticleType>
int EnergyRadws<ndim,ParticleType>::GetITemp
 (const FLOAT temp)
{
  int itemp = getClosestIndex(eos_temp, eos_temp + ntemp , temp);

  if (temp <= eos_temp[itemp]) {
    itemp = max(itemp - 1, 0);
  }

  return itemp;
}



//=================================================================================================
//  EnergyRadws::GetKappa()
/// GetKappa returns Kappa  for index of density and temp
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyRadws<ndim,ParticleType>::GetKappa
 (const int idens,
  const int itemp,
  const FLOAT logrho,
  const FLOAT logtemp,
  FLOAT &kappa,
  FLOAT &kappar,
  FLOAT &kappap)
{
	//cout << "getkappa1 " << idens << " " << itemp << " " << logrho << " " << logtemp << endl;
  const FLOAT epsilon = (logtemp - eos_temp[itemp])/(eos_temp[itemp+1]-eos_temp[itemp]);
  const FLOAT delta = (logrho - eos_dens[idens])/(eos_dens[idens+1]-eos_dens[idens]);

  kappa = kappa_table[idens+1][itemp+1]*epsilon*delta +
    kappa_table[idens+1][itemp]*(1-epsilon)*delta +
    kappa_table[idens][itemp+1]*epsilon*(1-delta) +
    kappa_table[idens][itemp]*(1-epsilon)*(1-delta);

  kappap = kappap_table[idens+1][itemp+1]*epsilon*delta +
    kappap_table[idens+1][itemp]*(1-epsilon)*delta +
    kappap_table[idens][itemp+1]*epsilon*(1-delta) +
    kappap_table[idens][itemp]*(1-epsilon)*(1-delta);

  kappar = kappap;

  //cout << "getkappa2 " << kappa << " " <<  kappap << endl;
  return;
}



//=================================================================================================
//  EnergyRadws::eBalance()
//  Calculates net heating rate due to hydro (i.e. expansion/contraction)
//  and radiative effects (i.e heating/cooling)
//=================================================================================================
template <int ndim, template <int> class ParticleType>
FLOAT EnergyRadws<ndim,ParticleType>::ebalance
 (const FLOAT dudt,
  const FLOAT temp_ex,
  const FLOAT temp,
  const FLOAT kappa,
  const FLOAT kappap,
  const FLOAT col2)
{
  return dudt - 4.0*rad_const*(pow(temp,4) - pow(temp_ex,4))/((col2*kappa) + (1.0/kappap));
}


//=================================================================================================
//  EnergyRadws::GetEnergy()
/// GetEnergy returns Energy  for index of density and temp
//=================================================================================================

template <int ndim, template <int> class ParticleType>
FLOAT EnergyRadws<ndim,ParticleType>:: GetEnergy(int idens, int itemp, FLOAT logrho, FLOAT logtemp)
{
  const FLOAT epsilon = (logtemp - eos_temp[itemp])/(eos_temp[itemp+1] - eos_temp[itemp]);
  const FLOAT delta = (logrho - eos_dens[idens])/(eos_dens[idens+1] - eos_dens[idens]);

  return eos_energy[idens+1][itemp+1]*epsilon*delta +
    eos_energy[idens+1][itemp]*(1 - epsilon)*delta +
    eos_energy[idens][itemp+1]*epsilon*(1 - delta) +
    eos_energy[idens][itemp]*(1 - epsilon)*(1 - delta);
 }



//=================================================================================================
//  EnergyRadws::GetMuBar()
/// GetMuBar returns MuBar  for index of density and temp
//=================================================================================================
template <int ndim, template <int> class ParticleType>
FLOAT EnergyRadws<ndim,ParticleType>::GetMuBar
 (const int idens,
  const int itemp,
  const FLOAT logrho,
  const FLOAT logtemp)
{
  const FLOAT epsilon = (logtemp - eos_temp[itemp])/(eos_temp[itemp+1] - eos_temp[itemp]);
  const FLOAT delta = (logrho - eos_dens[idens])/(eos_dens[idens+1] - eos_dens[idens]);

  return eos_mu[idens+1][itemp+1]*epsilon*delta +
    eos_mu[idens+1][itemp]*(1 - epsilon)*delta +
    eos_mu[idens][itemp+1]*epsilon*(1 - delta) +
    eos_mu[idens][itemp]*(1 - epsilon)*(1 - delta);
}



//=================================================================================================
//  EnergyRadws::GetTemp()
/// GetTemp returns the temperature of particle form it's density and specific internal energy.
/// The value is taken from the EoS table, not via a direct calculation.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
FLOAT EnergyRadws<ndim,ParticleType>::GetTemp
 (Particle<ndim> &part)
{
  FLOAT logrho = log10(part.rho);
  int idens = GetIDens(logrho);

  FLOAT result = max(part.temp_ambient,temp_ambient0);

  if (part.u == 0.0) return result;

  int temp_index = -1;
  FLOAT udif = 1E20;
  for (int i = 0; i < ntemp; ++i) {
    if (fabs(part.u - eos_energy[idens][i]) < udif) {
      udif = fabs(part.u - eos_energy[idens][i]);
      temp_index = i;
    }
  }

  if (temp_index > ntemp - 1) temp_index = ntemp - 1;

  result = pow10(eos_temp[temp_index]);

  //part.gamma = ((8.2673E3 * result * temp_unit) / (eos_mu[idens][temp_index] * eos_energy[idens][temp_index] * u_unit * 1E-4)) + 1.0;
  part.mu_bar = eos_mu[idens][temp_index];
  part.gamma = eos_gamma[idens][temp_index];
  //part.mu_bar = (FLOAT) 2.353;
  //part.gamma = (FLOAT) 5.0/3.0;
  return result;
}



template class EnergyRadws<1, GradhSphParticle>;
template class EnergyRadws<2, GradhSphParticle>;
template class EnergyRadws<3, GradhSphParticle>;
template class EnergyRadws<1, SM2012SphParticle>;
template class EnergyRadws<2, SM2012SphParticle>;
template class EnergyRadws<3, SM2012SphParticle>;
template class EnergyRadws<1, MeshlessFVParticle>;
template class EnergyRadws<2, MeshlessFVParticle>;
template class EnergyRadws<3, MeshlessFVParticle>;
