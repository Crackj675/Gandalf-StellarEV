//=================================================================================================
//  StellarEvolution.cpp
//  All routines for creating new sinks and accreting gas and updating all
//  sink particle propterties.
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
#include <limits>
#include "Precision.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "Parameters.h"
#include "Sph.h"
#include "Nbody.h"
#include "Sinks.h"
#include "Debug.h"
#include "StellarEvolution.h"
#include "Exception.h"
#include "InlineFuncs.h"
#include "Simulation.h"
using namespace std;





//=================================================================================================
//  StellarEvolution::StellarEvolution()
/// StellarEvolution class constructor
//=================================================================================================
template <int ndim>
StellarEvolution<ndim>::StellarEvolution(SimUnits *simunits_, SimulationBase *sim_,FLOAT AccFBminMass_, string acc_fb_)
{
  FLOAT num, denom, tempunit;
  simunits      = simunits_;
  sim           = sim_;
  acc_fb        = acc_fb_;

  if (acc_fb == "none"){
    cout << "Attetion: Use this model in combination with 'Accretion Feedback' and 'RadWS' for radiative Feedback!" << endl;
  }

  string line;
  FLOAT   tab_n_, tab_rad_, tab_DthetaDdx_, tab_rhoC_rhoBar_, tab_Bn_, tab_x2DthetaDdx_;
  FLOAT   tab_mass_, tab_beta_, tab_dlogbeta_, tab_dlogbbc_;

	cout << "Angmom Unit: " << simunits->angmom.outSI << endl;

  f_k   = 0.5;
  f_acc = 0.5;
  f_rad = 0.33;
	MinMass = AccFBminMass_;			///< Min. mass
  ntab1 = 100;
  ntab2 = 50;
  T_H   = 3000.0/(simunits->temp.outscale * simunits->temp.outSI);


  // compute rad_const = Stefan-Boltzmann in code units
  num       = pow(simunits->r.outscale*simunits->r.outSI,2)*simunits->t.outscale*simunits->t.outSI;
  denom     = simunits->E.outscale * simunits->E.outSI;
  tempunit  = simunits->temp.outscale * simunits->temp.outSI;
  rad_const = stefboltz*(num*pow(tempunit,4.0))/denom;

  // compute G_Const in code units
  num       = pow(simunits->t.outscale*simunits->t.outSI,2)*simunits->m.outscale*simunits->m.outSI;
  denom     = pow(simunits->r.outscale*simunits->r.outSI,3);
  G_Const   = G_const*num/denom;

  // compute L_sun in code units
  denom     = simunits->L.outscale*simunits->L.outSI;
  L_Sun     = L_sun/denom;

  // compute M_sun in code units
  denom     = simunits->m.outscale*simunits->m.outSI;
  M_Sun     = m_sun/denom;

  // compute R_sun in code units
  denom     = simunits->r.outscale*simunits->r.outSI;
  R_Sun     = r_sun/denom;

  // compute a_Const in code units
  num       = pow(simunits->temp.outscale*simunits->temp.outSI,4)*pow(simunits->r.outscale*simunits->r.outSI,3);
  denom     = simunits->E.outscale*simunits->E.outSI;
  a_Const   = 7.56591e-16*num/denom;
  //a_Const   = 2.5743845e-65*num/denom;
  //a_Const   = 1.673534e-37*num/denom;

  // compute k_B in code units
  num       = simunits->temp.outscale*simunits->temp.outSI;
  denom     = simunits->E.outscale*simunits->E.outSI;
  k_B       = k_boltzmann*num/denom;

  // compute m_H in code units
  denom     = simunits->m.outscale*simunits->m.outSI;
  m_H       = m_hydrogen/denom;


  // compute time unit
  denom     = simunits->t.outscale*simunits->t.outSI;
  yr0       = yr/denom;

  // read tables
  tab_n               = new FLOAT[ntab1];
  tab_rad             = new FLOAT[ntab1];
  tab_DthetaDdx       = new FLOAT[ntab1];
  tab_rhoC_rhoBar     = new FLOAT[ntab1];
  tab_Bn              = new FLOAT[ntab1];
  tab_x2DthetaDdx     = new FLOAT[ntab1];
  tab_mass            = new FLOAT[ntab1];
  tab_beta            = new FLOAT*[ntab2-2];
  tab_dlogbeta        = new FLOAT*[ntab2-2];
  tab_dlogbbc         = new FLOAT*[ntab2-2];

  for (int i=0; i<ntab2-2; i++){
      tab_beta[i]     = new FLOAT[ntab1];
      tab_dlogbeta[i] = new FLOAT[ntab1];
      tab_dlogbbc[i]  = new FLOAT[ntab1];

    }

  int i = 0;
  ifstream file;
  file.open("modeldata_table.dat", ios::in);
  while (getline(file, line)) {
    istringstream istr(line);

    if (istr >> tab_n_ >> tab_rad_ >> tab_DthetaDdx_ >> tab_rhoC_rhoBar_>> tab_Bn_ >> tab_x2DthetaDdx_ ) {
      tab_n[i]              = tab_n_;
      tab_rad[i]            = tab_rad_;
      tab_DthetaDdx[i]      = tab_DthetaDdx_;
      tab_rhoC_rhoBar[i]    = tab_rhoC_rhoBar_;
      tab_Bn[i]             = tab_Bn_;
      tab_x2DthetaDdx[i]    = tab_x2DthetaDdx_;
      ++i;
    }
  }
  file.close();

  i = 0;
  int j = 0;
  ifstream betafile;
  betafile.open("beta_table.dat", ios::in);
  while (getline(betafile, line)) {
    istringstream istr(line);
    if (istr >> tab_mass_ >> tab_beta_ >> tab_dlogbeta_ >> tab_dlogbbc_) {

      tab_mass[i]            = tab_mass_;
      tab_beta[i][j]         = tab_beta_;
      tab_dlogbeta[i][j]     = tab_dlogbeta_;
      tab_dlogbbc[i][j]      = tab_dlogbbc_;
      ++i;
      if (i == ntab2-2){
        i = 0;
        ++j;
      }
    }
  }
  betafile.close();

  return;
}

template <int ndim>
void StellarEvolution<ndim>::StellarEvolutionMain(Sinks<ndim> *sinks){

  // Stellar properties
  int stage;
  FLOAT dm;            // Accreted mass onto the star during this timestep      ToDo: check! + if no AccretionFB
  FLOAT dt;            // timestep                                              ToDo: check if right timestep
  FLOAT dr;            // Change of stellar radius
  FLOAT m;
  FLOAT n;
  FLOAT dmdt;
  FLOAT mdmdt;
  FLOAT r;
  FLOAT a_g;
  FLOAT mD;

  // Table data
  FLOAT rad;
  FLOAT DthetaDdx;
  FLOAT rhoC_rhoBar;
  FLOAT Bn;
  FLOAT beta;
  FLOAT x2DthetaDdx;
  FLOAT dlogbeta;
  FLOAT dlogbbc;

  // Luminosities
  FLOAT L_int;
  FLOAT L_ms;
  FLOAT L_I;
  FLOAT L_D;
  FLOAT L_H;
  FLOAT L_tot;

  Nsink = sinks->Nsink;

  for (int i=0; i<Nsink; i++) {

    // Update stellar properties
    if (acc_fb == "episodic"){
      dm    = sinks->sink[i].M_star - sinks->sink[i].M_star_old;
      m     = sinks->sink[i].M_star;
      //cout << "Stellar Mass : " << m << endl;
    } else {
      dm    = sinks->sink[i].star->m - sinks->sink[i].m_old;
      m     = sinks->sink[i].star->m;
    }
    dt    = sim->timestep*sinks->sink[i].star->nstep;
    dmdt  = dm/dt;
    mdmdt = (dmdt*yr0)/(M_Sun);
    mD    = sinks->sink[i].stellarMd;
    r     = sinks->sink[i].stellarRadius;
    n     = sinks->sink[i].stellarPI;
    stage = sinks->sink[i].stellarStage;
    a_g   = 3.0/(5.0 - n);
    //sinks->sink[i].dmdt = dmdt;

    //cout << "L_ZAMS " << ZAMSlum(1.0 / simunits->m.outscale)*simunits->L.outscale << endl;
    if(stage == 0){
      InitStar(&sinks->sink[i], m, mdmdt);
      mD    = sinks->sink[i].stellarMd;
      r     = sinks->sink[i].stellarRadius;
      n     = sinks->sink[i].stellarPI;
      stage = sinks->sink[i].stellarStage;
    }
    else if(stage == 5){
      r     = ZAMSradius(m);
      L_tot = ZAMSlum(m);
    }
    else if(stage > 0 and stage < 5){
      getModelData(n, rad, DthetaDdx, rhoC_rhoBar, Bn, x2DthetaDdx);
      getBetaData(n, m, beta, dlogbeta, dlogbbc);

      UpdateStellarEvolution(&sinks->sink[i], L_H, L_ms, L_int, L_I, L_D, r, mD,
        m, dm,  dmdt, mdmdt, dt, a_g, beta, dlogbeta, dlogbbc, stage);
      ComputeLumninosity(m, r, dmdt, L_int, L_tot);
      AdvanceEvolutionaryState(&sinks->sink[i], stage, r, n, m, mD, L_ms, L_D, rhoC_rhoBar, Bn);
      updateStellarProp(&sinks->sink[i], m, n, r, mD, L_tot, stage);
    }
	}
  return;
}

template <int ndim>
void StellarEvolution<ndim>::InitStar(SinkParticle<ndim> *sink, FLOAT m, FLOAT mdmdt){

  if (m > MinMass){

    //sink->stellarRadius   = 2.5*r_sun*pow(mdmdt/1e-5,0.2);
    sink->stellarRadius   = 2.5*R_Sun*pow(mdmdt/1e-5,0.2);
    sink->stellarPI       = min(max(5.0 - 3.0/(1.475 + 0.07*log10(mdmdt)),1.5),3.0);
    sink->stellarMd       = m;
    sink->stellarStage    = 1;
    //cout << "Radius : " << sink->stellarRadius/(r_sun) << endl;
    cout << "--------------------> Init stellear properties <--------------------" << endl;
    cout << "r " << sink->stellarRadius/R_Sun*simunits->r.outscale  << "   n " << sink->stellarPI << "   m_d " << sink->stellarMd*simunits->m.outscale << "   Stage " << sink->stellarStage << endl;
    cout << "--------------------------------------------------------------------" << endl;
  }
  return;
}

template <int ndim>
void StellarEvolution<ndim>::UpdateStellarEvolution(SinkParticle<ndim> *sink, FLOAT &L_H, FLOAT &L_ms, FLOAT &L_int,
  FLOAT &L_I, FLOAT &L_D, FLOAT &r, FLOAT &mD, FLOAT m, FLOAT dm, FLOAT dmdt, FLOAT mdmdt, FLOAT dt, FLOAT a_g, FLOAT beta, FLOAT dlogbeta, FLOAT dlogbbc, int stage){

  FLOAT dr;

  // Compute Luminosities
  L_H   = 4.0*pi*r*r*rad_const*pow(T_H,4);
  L_ms  = ZAMSlum(m);
  L_int = max(L_H, L_ms);
  L_I   = 2.5*1e5*L_Sun*mdmdt;
  computeLd(stage, L_int, L_I, m, r, dmdt, mdmdt, a_g, f_k, beta, dlogbbc, L_D);
  mD += dm - 1e-5*M_Sun*(L_D/(15*L_Sun)*(dt/yr0));

  //cout << " L_H " << L_H << " L_ms " << L_ms << " L_int "  << L_int << " L_I " << L_I << endl;
  // Update Radius
  dr = 2.0 * (dm/m) * (1.0 - (1.0-f_k)/(a_g*beta) + 0.5*dlogbeta)*r
      -2.0*dt/(a_g*beta)*(r/(G_Const*m*m))*(L_int + L_I - L_D)*r;
  r += dr;
  //cout << "Radius " << r*simunits->r.outscale *r_pc /r_sun << "     dr " << dr*simunits->r.outscale*r_pc /r_sun << endl;

  return;
}

template <int ndim>
void StellarEvolution<ndim>::ComputeLumninosity(FLOAT m, FLOAT r, FLOAT dmdt, FLOAT L_int, FLOAT &L_tot){

  FLOAT L_acc;
  FLOAT L_disk;

  L_acc = f_acc* f_k * G_Const * m * dmdt / r;
  L_disk = (1. - f_k) * G_Const * m * dmdt / r;
  L_tot = L_int + L_acc + L_disk;

  return;
}

template <int ndim>
void StellarEvolution<ndim>::AdvanceEvolutionaryState(SinkParticle<ndim> *sink, int &stage,
  FLOAT &r, FLOAT &n, FLOAT m, FLOAT mD, FLOAT L_ms, FLOAT L_D, FLOAT rhoC_rhoBar, FLOAT Bn){

  //cout << m << " " << mD << " " << L_D/L_ms << " " << stage << endl;
  switch(stage){
    case 1:
      FLOAT T_c;
      computeTc(m, r, Bn, rhoC_rhoBar, T_c);
      sink->stellarTc = T_c;
      if (T_c > 1.5e6/(simunits->temp.outscale*simunits->temp.outSI)){
        //cout << "T_c " << T_c << endl;
        cout << "STAGE 2 " << endl;
        n     = 1.5;
        stage  = 2;
      }
      break;
    case 2:
      if (mD <= 0.0){
        stage = 3;
        cout << "STAGE 3 " << endl;
      }
      break;
    case 3:
      //cout << "Lum " << L_D << " " << L_ms << " " << L_D/L_ms << endl;
      if (L_D/L_ms < f_rad){
        stage = 4;
        n     = 3.0;
        r    *= 2.1;

        cout << "STAGE 4 " << endl;
      }
      if (r <= ZAMSradius(m)){
        stage = 5;
        cout << "STAGE 5 " << endl;
      }
      break;
    case 4:
      if (r <= ZAMSradius(m)){
        stage = 5;
        cout << "STAGE 5 " << endl;
      }
      break;
  }
  return;
}

template <int ndim>
FLOAT StellarEvolution<ndim>::ZAMSradius(FLOAT m){

  FLOAT mass;
  FLOAT Num;
  FLOAT Denom;
  FLOAT a = 1.71535900;
  FLOAT b = 6.59778800;
  FLOAT c = 10.08855000;
  FLOAT d = 1.01249500;
  FLOAT e = 0.07490166;
  FLOAT f = 0.01077422;
  FLOAT g = 3.08223400;
  FLOAT h = 17.84778000;
  FLOAT i = 0.00022582;

  mass = m * simunits->m.outscale;
  Num     = a*pow(mass, 2.5) + b*pow(mass, 6.5) + c*pow(mass, 11) + d*pow(mass, 19) + e*pow(mass, 19.5);
  Denom   = f + g*pow(mass, 2) + h*pow(mass, 8.5) + pow(mass, 18.5) + i*pow(mass, 19.5);

  return r_sun * Num / (Denom * r_pc *simunits->r.outscale) ;
}

// ZAMS lumninosities fitted with Eq. (1) from Tout(1996)
template <int ndim>
FLOAT StellarEvolution<ndim>::ZAMSlum(FLOAT m){

  FLOAT mass;
  FLOAT Num;
  FLOAT Denom;
  FLOAT a = 0.39704170;
  FLOAT b = 8.52762600;
  FLOAT c = 0.00025546;
  FLOAT d = 5.43288900;
  FLOAT e = 5.56357900;
  FLOAT f = 0.78866060;
  FLOAT g = 0.00586685;
  mass = m * simunits->m.outscale;
  Num     = a*pow(mass, 5.5) + b*pow(mass, 11);
  Denom   = c + pow(mass, 3) + d*pow(mass, 5) + e*pow(mass, 7) + f*pow(mass, 8) + g*pow(mass, 9.5);

  return Num / (Denom * simunits->L.outscale);
}

template <int ndim>
void StellarEvolution<ndim>::getModelData(FLOAT n, FLOAT &rad, FLOAT &DthetaDdx,
  FLOAT &rhoC_rhoBar, FLOAT &Bn, FLOAT &x2DthetaDdx){

  int nIndex;
  FLOAT epsilon;

  if (n == 3) nIndex = ntab1-1;
  else{
    for (int i=0; i<ntab1-2; i++) {
      if (n >= tab_n[i] and n < tab_n[i+1]){
        nIndex = i;
        break;
      }
    }
  }

  epsilon     =  (n - tab_n[nIndex])/(tab_n[nIndex+1]-tab_n[nIndex]);

  rad         =  epsilon*tab_rad[nIndex]         + (1.0-epsilon)*tab_rad[nIndex+1];
  DthetaDdx   =  epsilon*tab_DthetaDdx[nIndex]   + (1.0-epsilon)*tab_DthetaDdx[nIndex+1];
  rhoC_rhoBar =  epsilon*tab_rhoC_rhoBar[nIndex] + (1.0-epsilon)*tab_rhoC_rhoBar[nIndex+1];
  Bn          =  epsilon*tab_Bn[nIndex]          + (1.0-epsilon)*tab_Bn[nIndex+1];
  x2DthetaDdx =  epsilon*tab_x2DthetaDdx[nIndex] + (1.0-epsilon)*tab_x2DthetaDdx[nIndex+1];
  return;
}

template <int ndim>
void StellarEvolution<ndim>::getBetaData(FLOAT n, FLOAT m, FLOAT &beta,
  FLOAT &dlogbeta, FLOAT &dlogbbc){

  int nIndex;
  int mIndex;
  FLOAT epsilon;
  FLOAT delta;


  if (n == 3) nIndex = ntab1-1;
  else{
    for (int i=0; i<ntab1-2; i++) {
      if (n >= tab_n[i] and n < tab_n[i+1]){
        nIndex = i;
        break;
      }
    }
  }
  for (int i=0; i<ntab2-4; i++) {
    if (m >= tab_mass[i] and m < tab_mass[i+1]){
      mIndex = i;
    }
  }

  epsilon     = (n - tab_n[nIndex])/(tab_n[nIndex+1]-tab_n[nIndex]);
  delta       = (m - tab_mass[mIndex])/(tab_mass[mIndex+1]-tab_mass[mIndex]);

  beta        = epsilon*delta*tab_beta[mIndex+1][nIndex+1] + (1.0-epsilon)*delta*tab_beta[mIndex+1][nIndex] + epsilon*(1.0-delta)*tab_beta[mIndex][nIndex+1] + (1.0-epsilon)*(1.0-delta)*tab_beta[mIndex][nIndex];
  dlogbeta    = epsilon*delta*tab_dlogbeta[mIndex+1][nIndex+1] + (1.0-epsilon)*delta*tab_dlogbeta[mIndex+1][nIndex] + epsilon*(1.0-delta)*tab_dlogbeta[mIndex][nIndex+1] + (1.0-epsilon)*(1.0-delta)*tab_dlogbeta[mIndex][nIndex];
  dlogbbc     = epsilon*delta*tab_dlogbbc[mIndex+1][nIndex+1] + (1.0-epsilon)*delta*tab_dlogbbc[mIndex+1][nIndex] + epsilon*(1.0-delta)*tab_dlogbbc[mIndex][nIndex+1] + (1.0-epsilon)*(1.0-delta)*tab_dlogbbc[mIndex][nIndex];

  //cout << "mIndex " << mIndex << " nIndex " << nIndex << " beta " << beta << " dlogbeta "  << dlogbeta << " dlogbbc " << dlogbbc << endl;

  return;
}

template <int ndim>
void StellarEvolution<ndim>::computeLd(int stage, FLOAT L_int, FLOAT L_I,
  FLOAT m, FLOAT r, FLOAT dmdt, FLOAT mdmdt, FLOAT a_g, FLOAT f_k, FLOAT beta, FLOAT dlogbbc, FLOAT &L_D){

  if (stage < 2){
    L_D = 0.;
  }
  else if (stage == 2){
    L_D = L_int + L_I + (G_Const * m)/r * dmdt * ( 1. - f_k - a_g*beta/2. * (1. + dlogbbc));
      //cout << "L_D " << L_D << " " << stage <<  " " << dmdt << endl;
  }
  else if (stage > 2){
    L_D = 15. * L_Sun * mdmdt * 1.e5;
      //cout << "L_D " << L_D << " " << stage <<  " " << dmdt << endl;
  }
  return;
}

template <int ndim>
void StellarEvolution<ndim>::computeTc(FLOAT m, FLOAT r, FLOAT Bn, FLOAT rhoC_rhoBar, FLOAT &T_c){
// Solving the quartic using Ferrari's method.


  FLOAT rho_avg;    //rho average
  FLOAT rho_c;      //rho center
  FLOAT e;
  FLOAT a;
  FLOAT d;
  FLOAT q;
  FLOAT delta0;
  FLOAT delta1;
  FLOAT Q;
  FLOAT S;

  rho_avg = 3. * m / (4. * pi * pow(r,3));
  rho_c   = rhoC_rhoBar * rho_avg;
  e       = -pow((4. * pi),(1./3.)) * Bn * G_Const * pow(m,(2./3.)) * pow(rho_c,(4./3.));
  a       = (1./3.)*a_Const;
  d       = (rho_c*k_B)/(0.613*m_H);
  q       = d/a;
  delta0  = 12.0*a*e;
  delta1  = 27.0*a*d*d;
  //Q       = pow(0.5*(delta1+sqrt(27*(256*pow(a,3)*pow(e,3) - 27*a*a+pow(d,4)))),1./3.);
  Q       = pow(0.5*(delta1+sqrt(delta1*delta1-4.0*pow(delta0,3))),1./3.);
  S       = 0.5*sqrt(1./(3.*a)*(Q+delta0/Q));

  T_c = -S + 0.5*sqrt(-4.*S*S+q/S);

  //cout << (-S + 0.5*sqrt(-4.*S*S+q/S))*simunits->temp.outscale << endl;
  //cout << (-S - 0.5*sqrt(-4.*S*S+q/S))*simunits->temp.outscale << endl;
  //cout << ( S + 0.5*sqrt(-4.*S*S-q/S))*simunits->temp.outscale << endl;
  //cout << ( S - 0.5*sqrt(-4.*S*S-q/S))*simunits->temp.outscale << endl;

  //cout << "Tc : " << T_c*simunits->temp.outscale << endl;
  return;
}

template <int ndim>
FLOAT StellarEvolution<ndim>::corepress(FLOAT Pc, FLOAT rho_c, FLOAT Tc){

  FLOAT corepress;
  corepress = (rho_c * k_B * Tc)/(0.613*m_H) + (1./3.)*a_Const*pow(Tc,4) - Pc;

  return corepress;
}
/*
template <int ndim>
void StellarEvolution<ndim>::computeTc(FLOAT m, FLOAT r, FLOAT Bn, FLOAT rhoC_rhoBar, FLOAT &T_c){

  FLOAT     T1,T2,tol,Pc;
  FLOAT     rt;
  FLOAT     func;
  int       maxiter=100;
  FLOAT     eps=std::numeric_limits<FLOAT>::epsilon();
  int       iter;
  FLOAT     a,b,c,d,e,fa,fb,fc,p,q,R,s,tol1,xm;
  FLOAT     rho_avg, rho_c;

  T1      = 0.0;
  T2      = 1e8/simunits->temp.outscale;
  tol     = 0.001;
  a       = T1;
  b       = T2;

  iter=0;

  rho_avg = 3. * m / (4. * pi * pow(r,3));
  rho_c   = rhoC_rhoBar * rho_avg;
  Pc      = pow((4. * pi),(1./3.)) * Bn * G_Const * pow(m,2./3.) * pow(rho_c,4./3.);

  c=a;
  fc=fa;
  d=b-a;
  e=d;

  fa = corepress(Pc,rho_c,a);
  fb = corepress(Pc,rho_c,b);

  for (iter=0; iter<maxiter; iter++) {
    iter=iter+1;

    if (fb*fc>0){
        c=a; fc=fa; d=b-a; e=d;
    }

    if (abs(fc)<abs(fb)){
        a=b; b=c; c=a;
        fa=fb; fb=fc; fc=fa;
    }

    tol=2*eps*abs(b)+t; m=(c-b)/2; //Toleranz

    if ((abs(m)>tol) and (abs(fb)>0)){ //Verfahren muss noch durchgeführt werden

      if ((abs(e)<tol) or (abs(fa)<=abs(fb))){
        d=m;
        e=m;
      }
      else{
        s=fb/fa;
        if (a==c){
          p=2*m*s;
          q=1-s;
        }
        else{
          q=fa/fc;
          R=fb/fc;
          p=s*(2*m*q*(q-R)-(b-a)*(R-1));
          q=(q-1)*(R-1)*(s-1);
        }
        if (p>0) q=-q;
        else{
          p=-p;
        }
        s=e;
        e=d;
        if (( 2*p<3*m*q-abs(tol*q) ) and (p<abs(s*q/2))){
          d=p/q;
        }
        else{
          d=m;
          e=m;
        }
      }

      a=b;
      fa=fb;

      if (abs(d)>tol){
          b=b+d;
      }
      else{
        if (m>0) b=b+tol;
        else b=b-tol;
      }
    }
    else{
      break;
    }
    fb=corepress(Pc,rho_c,b);
  }

  T_c = b;
  cout << "T_c : " << T_c *simunits->temp.outscale << endl;
  return;
}
*/
template <int ndim>
void StellarEvolution<ndim>::updateStellarProp(SinkParticle<ndim> *sink, FLOAT m, FLOAT n, FLOAT r, FLOAT mD, FLOAT L_tot, int stage){

  sink->stellarRadius   = r;
  sink->stellarPI       = n;
  sink->stellarMd       = mD;
  sink->stellarStage    = stage;
  sink->luminosity      = L_tot;
  return;
}

template class StellarEvolution<1>;
template class StellarEvolution<2>;
template class StellarEvolution<3>;










