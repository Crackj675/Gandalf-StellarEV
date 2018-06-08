//=================================================================================================
//  EnergyEquation.h
//  Class definitions of main energy equation class plus inherited children
//  classes for various energy integration algorithms.
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
#ifndef _ACCRETIONFB_H_
#define _ACCRETIONFB_H_

#include "Constants.h"
#include "Particle.h"
#include "Precision.h"
#include "SimUnits.h"
#include "Hydrodynamics.h"
#include "Sinks.h"
#include "Nbody.h"

//=================================================================================================
//  AccretionFB
/// Main energy equation class, with virtual functions that require full
/// definitions in the children classes.
//=================================================================================================
template <int ndim>
class AccretionFB
{
 public:



  	AccretionFB(SimUnits *, SimulationBase *, string, const FLOAT, const FLOAT, const FLOAT);
  	~AccretionFB();

  	virtual void Luminosity(Hydrodynamics<ndim> *, Sinks<ndim> *, const FLOAT)= 0;
  	void updateAccFBTemp(const int, Hydrodynamics<ndim> *, Sinks<ndim> *);
  	virtual void Accretion(Sinks<ndim> *, Hydrodynamics<ndim> *) = 0;

  string stellarEv;
	FLOAT rad_const;
	FLOAT L_Sun;
	FLOAT M_Sun;
	FLOAT G_Const;
	FLOAT R_Sun;
	FLOAT mu_bar;
	FLOAT gammam1;
	FLOAT temp_ambient0;

  	Sinks<ndim> *sinks;
  	SimUnits *simunits;
    SimulationBase *sim;

  	//Hydrodynamics<ndim> *hydro;

  //SinkParticle<ndim> *sinkdata;
  //Hydrodynamics<ndim> *hydro;          ///< Hydrodynamics algorithm pointer
  //Nbody<ndim> *nbody;
};

//=================================================================================================
//  ContinuousAFB
/// Continuous Accretion Feedback
//=================================================================================================
template <int ndim>
class ContinuousAFB : public AccretionFB<ndim>
{

  using AccretionFB<ndim>::rad_const;
  using AccretionFB<ndim>::L_Sun;
  using AccretionFB<ndim>::M_Sun;
  using AccretionFB<ndim>::G_Const;
  using AccretionFB<ndim>::R_Sun;
  using AccretionFB<ndim>::simunits;
  using AccretionFB<ndim>::sinks;
  using AccretionFB<ndim>::sim;
  using AccretionFB<ndim>::stellarEv;

 public:

  ContinuousAFB(SimUnits *, SimulationBase *, string, const FLOAT, const FLOAT, const FLOAT, 
                const FLOAT, const FLOAT, const FLOAT);
  ~ContinuousAFB();

  void Luminosity(Hydrodynamics<ndim> *, Sinks<ndim> *, const FLOAT);
  void Accretion(Sinks<ndim> *, Hydrodynamics<ndim> *);

  FLOAT AccFBminMass;
  FLOAT f_angmom;
  FLOAT dmdt_max;
};

//=================================================================================================
//  EpisodicAFB
/// Episodic Accretion Feedback
//=================================================================================================
template <int ndim>
class EpisodicAFB : public AccretionFB<ndim>
{

  using AccretionFB<ndim>::rad_const;
  using AccretionFB<ndim>::L_Sun;
  using AccretionFB<ndim>::M_Sun;
  using AccretionFB<ndim>::G_Const;
  using AccretionFB<ndim>::R_Sun;
  using AccretionFB<ndim>::simunits;
  using AccretionFB<ndim>::sinks;
  using AccretionFB<ndim>::sim;
  using AccretionFB<ndim>::stellarEv;

 public:

  EpisodicAFB(SimUnits *, SimulationBase *, string, const FLOAT, const FLOAT, const FLOAT,
  				const FLOAT, const FLOAT, const FLOAT, const FLOAT);
  ~EpisodicAFB();

  void Luminosity(Hydrodynamics<ndim> *, Sinks<ndim> *, const FLOAT);
  void Accretion(Sinks<ndim> *, Hydrodynamics<ndim> *);

  Hydrodynamics<ndim> *hydro;
  float AccFBminMass;
  FLOAT alpha_MRI;
  FLOAT dmdt_BRG;
  FLOAT yr;
  FLOAT dmdt_MRI;
  FLOAT dm;
  FLOAT f_angmom;
};


#endif
