//=================================================================================================
//  Outflow.h
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
#ifndef _OUTFLOWFB_H_
#define _OUTFLOWFB_H_

#include "Constants.h"
#include "Particle.h"
#include "Precision.h"
#include "SimUnits.h"
#include "Hydrodynamics.h"
#include "Sinks.h"
#include "Nbody.h"
//#include "AccretionFB.h"

//=================================================================================================
//  Outflow
/// ...
/// ...
//=================================================================================================
template <int ndim>
class Outflow
{
 public:

  Outflow(SimUnits *,SimulationBase *, RandomNumber *, string, int, const FLOAT, const FLOAT,
          const FLOAT, const FLOAT, const FLOAT, const FLOAT, const FLOAT);
  ~Outflow();

	void Check_Outflow(Sinks<ndim> *, Hydrodynamics<ndim> *, const FLOAT, const int, const int, const int);
	void Calc_Outflow(Sinks<ndim> *, Hydrodynamics<ndim> *, const FLOAT, const int, const int, const int,
          const int, const int, const FLOAT);
	void Set_PV(int, int, int, FLOAT, FLOAT *, FLOAT *, FLOAT *);
	void Set_part(int, int, int, int, int, FLOAT, FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *, Hydrodynamics<ndim> *);
	void Update_sink(int, int, FLOAT*);
	void CallOutput(int, int, FLOAT*, FLOAT*, FLOAT);
  FLOAT ComputeAngmom(int, int);
	FLOAT dot_product(FLOAT* , FLOAT* );
	FLOAT* cross_product(FLOAT* , FLOAT*, FLOAT*);
	FLOAT* QuaternionMirror(FLOAT*, FLOAT*);
	FLOAT* Hprod(FLOAT* , FLOAT*, FLOAT* );
	FLOAT norm(FLOAT *);

	SimulationBase *sim;
	RandomNumber *randnumb;
	SimUnits *simunits;
	Sinks<ndim> *sinks;


	bool outflowTiming;
	int Nout;
  string stellarEv;
	FLOAT G_Const;
	FLOAT op_angle;
	FLOAT rad_const;
	FLOAT f_vel;
	FLOAT R_launch;
 	FLOAT theta0;
	FLOAT f_angmom;
	FLOAT f_eject;
	FLOAT r_eject;
	FLOAT dm;

};
#endif
