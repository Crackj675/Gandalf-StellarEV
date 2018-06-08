//=================================================================================================
//  StellarEvolution.h
//  Main stellar evolution class
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


#ifndef _StellarEvolution_H_
#define _StellarEvolution_H_


#include <string>
#include "Precision.h"
#include "CodeTiming.h"
#include "Constants.h"
#include "Hydrodynamics.h"
#include "Sinks.h"
#include "Nbody.h"
#include "NbodyParticle.h"
#include "NeighbourSearch.h"
#include "Parameters.h"
#include "Particle.h"
#include "SmoothingKernel.h"
#include "StarParticle.h"




//=================================================================================================
//  Class StellarEvolution
/// \brief   Main stellar evolution class.
/// \details ...
/// \author  P. Rohde
/// \date    21/03/2018
//=================================================================================================
template<int ndim>
class StellarEvolution
{
 public:

	// Constructor and destructor
	//-----------------------------------------------------------------------------------------------
  StellarEvolution(SimUnits *, SimulationBase *, FLOAT, string);
	~StellarEvolution();


	// Function prototypes
	//-----------------------------------------------------------------------------------------------
	void StellarEvolutionMain(Sinks<ndim> *);
  void InitStar(SinkParticle<ndim> *, FLOAT, FLOAT);
	void UpdateStellarEvolution(SinkParticle<ndim> *, FLOAT &, FLOAT &, FLOAT &, FLOAT &, FLOAT &, FLOAT &,
    FLOAT &, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, int);
	void ComputeLumninosity(FLOAT, FLOAT, FLOAT, FLOAT, FLOAT &);
	void AdvanceEvolutionaryState(SinkParticle<ndim> *, int &, FLOAT &, FLOAT &, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT);
  void getModelData(FLOAT, FLOAT &, FLOAT &, FLOAT &, FLOAT &, FLOAT &);
  void getBetaData(FLOAT, FLOAT, FLOAT &, FLOAT &, FLOAT &);
  void computeLd(int, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT &);
  void computeTc(FLOAT, FLOAT, FLOAT, FLOAT, FLOAT &);
  void updateStellarProp(SinkParticle<ndim> *, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, int);
  FLOAT ZAMSradius(FLOAT);
  FLOAT ZAMSlum(FLOAT);
  FLOAT corepress(FLOAT, FLOAT, FLOAT);

  // Local pointers
  //-----------------------------------------------------------------------------------------------
  SimUnits *simunits;
  SimulationBase *sim;

	// Local class variables
	//-----------------------------------------------------------------------------------------------
	int Nsink; 					///< # Sink particles
  int ntab1;
  int ntab2;
  string acc_fb;
	FLOAT MinMass; 				///< Min. mass
  FLOAT rad_const;
  FLOAT L_Sun;
  FLOAT M_Sun;
  FLOAT G_Const;
  FLOAT a_Const;
  FLOAT R_Sun;
  FLOAT k_B;
  FLOAT m_H;
  FLOAT yr0;
  FLOAT f_k;
  FLOAT f_acc;
  FLOAT f_rad;
  FLOAT T_H;
  FLOAT t;

  FLOAT *tab_n;
  FLOAT *tab_rad;
  FLOAT *tab_DthetaDdx;
  FLOAT *tab_rhoC_rhoBar;
  FLOAT *tab_Bn;
  FLOAT *tab_x2DthetaDdx;
  FLOAT *tab_mass;
  FLOAT **tab_beta;
  FLOAT **tab_dlogbeta;
  FLOAT **tab_dlogbbc;
};
#endif
