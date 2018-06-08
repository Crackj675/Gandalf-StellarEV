//===========================================
//======================================================
//  Outflow.cpp
//  Contains functions for Outflow class.
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
#include "RandomNumber.h"
#include "Precision.h"
#include "Constants.h"
#include "Debug.h"
#include "EnergyEquation.h"
#include "Hydrodynamics.h"
#include "Parameters.h"
#include "Particle.h"
#include "SimUnits.h"
#include "Simulation.h"
#include "Outflow.h"
using namespace std;



//=================================================================================================
//  Outflows::Outflow()
/// Outflow class constructor
//=================================================================================================
template <int ndim>
Outflow<ndim>::Outflow(
  SimUnits *simunits_,
  SimulationBase *sim_,
  RandomNumber *randnumb_,
  string stellarEv_,
  const int Nout_,
	const FLOAT R_launch_,
  const FLOAT f_vel_,
  const FLOAT op_angle_,
  const FLOAT theta0_,
	const FLOAT f_angmom_,
  const FLOAT f_eject_,
  const FLOAT r_eject_)
{
	Nout = Nout_;
	FLOAT denom;
	FLOAT num;

	sim 				= sim_;
	randnumb 		= randnumb_;
	simunits 		= simunits_;
	op_angle 		= op_angle_;
  f_angmom 		= f_angmom_;
  f_eject 		= f_eject_;
	r_eject			= r_eject_;
	f_vel 			= f_vel_;
	theta0 			= theta0_;
  stellarEv   = stellarEv_;


	cout << " -------------------------------------------------" << endl;
	cout << "You are using the Outflow-module!!" << endl;
	cout << " -------------------------------------------------" << endl;
	if (ndim < 3){
		cout << " -------------------------------------------------" << endl;
		cout << "Outflow-model only works in 3D!" << endl;
		cout << " -------------------------------------------------" << endl;
	}

//-------------------------------------------------------------------------------------------
// Unit scaling
//-------------------------------------------------------------------------------------------

	//Ensure that outflows are started at nresync, but can continue afterwars
	outflowTiming = false;

	// compute radius of the protostar in code units
	denom 			= simunits->r.outscale*simunits->r.outSI;
  R_launch 			= R_launch_ * r_sun/denom;
	cout << "R_launch: " << R_launch*simunits->r.outscale << endl;

  // compute G_Const in code units
  num       		= pow(simunits->t.outscale*simunits->t.outSI,2)*simunits->m.outscale*simunits->m.outSI;
  denom     		= pow(simunits->r.outscale*simunits->r.outSI,3);
  G_Const			= G_const*num/denom;


	r_eject /= 206000*simunits->r.outscale;
}

//=================================================================================================
//  Outflows::Outflow()
/// Outflow class destructor
//=================================================================================================
template <int ndim>
Outflow<ndim>::~Outflow()
{
}

//=================================================================================================
//  Outflows::Check_Outflow()
/// Main Outflow routine
//=================================================================================================
template <int ndim>
void Outflow<ndim>::Check_Outflow(Sinks<ndim> *sinks_, Hydrodynamics<ndim> *hydro_,
			const FLOAT current_time, const int n, const int level_step, const int level_max)
{
	sinks = sinks_;								// Sink particle pointer
	int 			Nsink = sinks->Nsink; 		// Number of sink-particle
	int 	 		i;


#pragma omp parallel for default(none) private(i) shared(Nsink, hydro_, cout)
	for (i=0; i<Nsink; i++) {

		int 		Npart;						// 1/4 of the amount of Particle for the outflow
		FLOAT 	dm_of; 						// Delta Mass for the outflow

		if (sinks->sink[i].acc_flag == 1){
			Particle<ndim>& part = hydro_->GetParticlePointer(1);
			dm = part.m;
			dm_of = f_eject*(sinks->sink[i].M_star - sinks->sink[i].M_star_old) + sinks->sink[i].dm_store;
			cout << dm_of << endl;
			if (dm_of > 4*Nout*dm){
				Npart = floor( 0.25 * dm_of / dm);
				dm_of -= 4*Npart*dm;
				Calc_Outflow(sinks, hydro_, current_time, Npart, i, n, level_step, level_max, dm);
			}
			else if (dm_of > 40*dm) {
				cout << "Attetion! --- At least 40 partciles are stored in the IAD! "<< endl;
			}
		sinks->sink[i].dm_store = dm_of;
		}
	}
	return;
}

//=================================================================================================
//  Outflows::Calc_Outflow()
/// Main Outflow routine
//=================================================================================================
template <int ndim>
void Outflow<ndim>::Calc_Outflow(Sinks<ndim> *sinks_, Hydrodynamics<ndim> *hydro_, const FLOAT current_time,
	const int Npart, const int SinkID, const int n, const int level_step, const int level_max, const FLOAT dm)
{
	int 	 			i,j;
	FLOAT				v0;               // Kepler velocity at R_out
	FLOAT 			L_norm;						// Angular Momentum of the IAD
	FLOAT 			L_total[ndim];    // Array to store total ang. mom.
	FLOAT 			e_L[ndim];        // normalized ang. mom. vec. of the IAD


	FLOAT 			pos[6*Npart];				// Position array
	FLOAT 			vel[6*Npart];				// Velovity array
	FLOAT 			pos_rev[6*Npart];   // Mirrored position array
	FLOAT 			vel_rev[6*Npart];   // Mirrored velocity array

	sinks = sinks_;								  // Sink particle pointer

//-------------------------------------------------------------------------------------------
// Read in information from sink particle
//-------------------------------------------------------------------------------------------

	//Update inner accretion Disc (IAD)
	for (int k=0; k<ndim; k++) sinks->sink[SinkID].L_IAD[k] = sinks->sink[SinkID].angmom[k] - sinks->sink[SinkID].L_star[k];
	for (j=0; j<ndim; j++) L_total[j] = 0.0;

	L_norm = norm(sinks->sink[SinkID].L_IAD);
	for (j=0; j<ndim; j++) e_L[j] = sinks->sink[SinkID].L_IAD[j] / L_norm;


// Calculate velocity scale
//-------------------------------------------------------------------------------------------

  if (stellarEv == "none"){
    v0 = sqrt((G_Const*sinks->sink[SinkID].star->m)/(R_launch));
  }
  else{
    v0 = sqrt((G_Const*sinks->sink[SinkID].star->m)/(sinks->sink[SinkID].stellarRadius));
  }
	if (R_launch > sinks->sink[SinkID].radius || sinks->sink[SinkID].stellarRadius > sinks->sink[SinkID].radius){
		v0 = sqrt((G_Const*sinks->sink[SinkID].M_star)/(sinks->sink[SinkID].radius));
		cout << "R_launch > R_sink : " << R_launch << " " << sinks->sink[SinkID].radius << endl;
	}

//-------------------------------------------------------------------------------------------
// Call function to set velocity and position for each outflow particle
//-------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(j) shared(e_L, v0, pos, vel)
	for (j=0; j<Npart; j++){
		Set_PV(j, Npart, SinkID, v0, pos, vel, e_L);
		assert(!(isnan(vel[j])));
	}

//-------------------------------------------------------------------------------------------
// Add symmetric particle
//-------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(j) shared(vel, pos, vel_rev, pos_rev)
	for (j=0; j<2*Npart; j++){
		vel_rev[ndim*j+0] = -vel[ndim*j+0];
		vel_rev[ndim*j+1] = -vel[ndim*j+1];
		vel_rev[ndim*j+2] = -vel[ndim*j+2];
		pos_rev[ndim*j+0] = -pos[ndim*j+0];
		pos_rev[ndim*j+1] = -pos[ndim*j+1];
		pos_rev[ndim*j+2] = -pos[ndim*j+2];
	}

	L_total[0] = 0.0;
	L_total[1] = 0.0;
	L_total[2] = 0.0;

	for (j=0; j<2*Npart; j++){
		L_total[0] += dm*(pos[ndim*j+1]*vel[ndim*j+2] - pos[ndim*j+2]*vel[ndim*j+1]);
		L_total[1] += dm*(pos[ndim*j+2]*vel[ndim*j+0] - pos[ndim*j+0]*vel[ndim*j+2]);
		L_total[2] += dm*(pos[ndim*j+0]*vel[ndim*j+1] - pos[ndim*j+1]*vel[ndim*j+0]);
		L_total[0] += dm*(pos_rev[ndim*j+1]*vel_rev[ndim*j+2] - pos_rev[ndim*j+2]*vel_rev[ndim*j+1]);
		L_total[1] += dm*(pos_rev[ndim*j+2]*vel_rev[ndim*j+0] - pos_rev[ndim*j+0]*vel_rev[ndim*j+2]);
		L_total[2] += dm*(pos_rev[ndim*j+0]*vel_rev[ndim*j+1] - pos_rev[ndim*j+1]*vel_rev[ndim*j+0]);
	}


//-------------------------------------------------------------------------------------------
// Transform into sink-particle reference system
//-------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(j) shared(vel, pos, vel_rev, pos_rev)
	for (j=0; j<2*Npart; j++){
		vel[ndim*j+0] += sinks->sink[SinkID].star->v[0];
		vel[ndim*j+1] += sinks->sink[SinkID].star->v[1];
		vel[ndim*j+2] += sinks->sink[SinkID].star->v[2];
		pos[ndim*j+0] += sinks->sink[SinkID].star->r[0];
		pos[ndim*j+1] += sinks->sink[SinkID].star->r[1];
		pos[ndim*j+2] += sinks->sink[SinkID].star->r[2];
		vel_rev[ndim*j+0] += sinks->sink[SinkID].star->v[0];
		vel_rev[ndim*j+1] += sinks->sink[SinkID].star->v[1];
		vel_rev[ndim*j+2] += sinks->sink[SinkID].star->v[2];
		pos_rev[ndim*j+0] += sinks->sink[SinkID].star->r[0];
		pos_rev[ndim*j+1] += sinks->sink[SinkID].star->r[1];
		pos_rev[ndim*j+2] += sinks->sink[SinkID].star->r[2];

		assert(!(isnan(vel[ndim*j+0])) ||!(isnan(vel[ndim*j+1])) || !(isnan(vel[ndim*j+2])));
		assert(!(isnan(pos[ndim*j+0])) ||!(isnan(pos[ndim*j+1])) || !(isnan(pos[ndim*j+2])));
	}

	cout << " ------------------------------------------------------------------" << endl;
	Set_part(SinkID, Npart, n, level_step, level_max, v0, current_time, pos, vel, pos_rev ,vel_rev, hydro_);
	Update_sink(SinkID, Npart, L_total);
	CallOutput(SinkID, Npart, L_total, e_L, v0);
	cout << " ------------------------------------------------------------------" << endl;

	return;
}


//=================================================================================================
//  Outflows::Set_PV()
/// Draws particle in within the ouflow-cone, transferres them to the angular momentum axis
/// and calculates the velocity components for the main and first mirrored particle
//=================================================================================================
template <int ndim>
void Outflow<ndim>::Set_PV(
	int part_ID,
	int Npart,
	int sink_ID,
	FLOAT V_part,
	FLOAT *pos,
	FLOAT *vel,
	FLOAT *e_L)
{

	int       j;
	FLOAT 		phi;				// azimuth angle
	FLOAT 	 	theta;				// polar angle
	FLOAT		  costheta;			// cos(theta)
  FLOAT 		costheta_min;
	FLOAT 		sintheta;			// sin(theta)
	FLOAT 		r_norm;				// distance

	FLOAT		  rotfreq;			// Angular frequency
	FLOAT 		f_vel_dist;			// velocity modulation
	FLOAT 		R_perp;				// perpendicular distance
	FLOAT		  V_rot_norm;			// norm of rot. vel. comp.
	FLOAT		  L_norm;				// total ang. mom. of both part.
	FLOAT 		z2;
  FLOAT     dL;

	FLOAT 		r[ndim];			// position array
	FLOAT 		r_sym[ndim];		// mirrored position array
	FLOAT 		V_rot[ndim];		// velocity array
	FLOAT 		V_rot_sym[ndim];	// mirrored velocity array
	FLOAT 		L_check[ndim];		// needed to sum up the ang. mom.
	FLOAT 		L_check_sym[ndim];	// needed to sum up the ang. mom.
	FLOAT 		L_check2[ndim];		// delete again!

//-------------------------------------------------------------------------------------------
// Draw phi and thetha within the opening angle
//-------------------------------------------------------------------------------------------

	// Draw costheta coressponing to the distribtion:
	// P(x) = sqrt(1.0/(1.0+theta0**2-x**2))/arctan(1./theta0)
	// (Matzner, 1999) with costheta < cos(op_angle)
	theta = op_angle;
	while (theta >= op_angle){
		z2 = pow(tan(randnumb->floatrand()*atan(1.0/theta0)),2);
		costheta = sqrt(z2*(theta0*theta0+1.0)/(1.0+z2));
		theta = acos(costheta);
	}
	sintheta = sin(theta);
	//r_norm	 = 2*sinks->sink[sink_ID].radius*randnumb->floatrand() + r_eject*sinks->sink[sink_ID].radius;
  r_norm   = (2.0*r_au / (r_pc*simunits->r.outscale))*sinks->sink[sink_ID].radius*randnumb->floatrand() + r_eject;
	phi      = (FLOAT) 2.0*pi*randnumb->floatrand();

	cout << "costheta " << costheta << " theta " << theta << " " << theta0 <<  endl;

// Draw phi and thetha within the opening angle
//-------------------------------------------------------------------------------------------
	r[0] 				= r_norm*sintheta*cos(phi);
	r[1] 				= r_norm*sintheta*sin(phi);
	r[2] 				= r_norm*costheta;
	r_sym[0] 		= r_norm*sintheta*cos(phi+pi);
	r_sym[1] 		= r_norm*sintheta*sin(phi+pi);
	r_sym[2] 		= r_norm*costheta;

// Compute rotational velocity component
//-------------------------------------------------------------------------------------------
  dL = ComputeAngmom(sink_ID, Npart);
	R_perp 				= r_norm *sintheta;
	rotfreq 			= dL/(dm*R_perp*R_perp);

	V_rot[0] 			= -rotfreq*r[1];
	V_rot[1] 			= rotfreq*r[0];
	V_rot[2] 			= 0.0;
	V_rot_sym[0] 	= -rotfreq*r_sym[1];
	V_rot_sym[1] 	= rotfreq*r_sym[0];
	V_rot_sym[2] 	= 0.0;


// Threshold for the rotational velocity component
//-------------------------------------------------------------------------------------------
	V_rot_norm 		= norm(V_rot);
	if (5*V_rot_norm > V_part){
		V_rot[0] 		*= V_part/(5*V_rot_norm);
		V_rot[1] 		*= V_part/(5*V_rot_norm);
		V_rot_sym[0] 	*= V_part/(5*V_rot_norm);
		V_rot_sym[1] 	*= V_part/(5*V_rot_norm);
		cout << " Rotational > radial velocity component:  " << V_part/(5*V_rot_norm) << endl;
	}


//-------------------------------------------------------------------------------------------
// Calculate position and rotational velocity component by rotating to angular momentum axis
//-------------------------------------------------------------------------------------------

	QuaternionMirror(e_L, r);
	QuaternionMirror(e_L, r_sym);
	QuaternionMirror(e_L, V_rot);
	QuaternionMirror(e_L, V_rot_sym);

	for (j=0; j<ndim; j++){
		pos[3 *  part_ID + j] 			= r[j];
		pos[3 * (part_ID+Npart) + j] 	= r_sym[j];
	}


//-------------------------------------------------------------------------------------------
// Calculate total velocity of both particle
//-------------------------------------------------------------------------------------------
	f_vel_dist = f_vel*sqrt(1.0/(atan(1.0/theta0)* sqrt(1.0+theta0*theta0-costheta*costheta)));
  for (j=0; j<ndim; j++){
		vel[3 *  part_ID + j]	 		=  f_vel_dist *r[j]* 	  V_part / r_norm; //f_vel_dist *
		vel[3 * (part_ID+Npart) + j] 	=  f_vel_dist *r_sym[j]* V_part / r_norm; //f_vel_dist *
	}


  //ofstream Outflow;
  //Outflow.open ("./Outflow.txt", std::ios_base::app);
  //Outflow << pos[3 *  part_ID + 0]*simunits->r.outscale << " " << pos[3 *  part_ID + 1]*simunits->r.outscale << " " << pos[3 *  part_ID + 2]*simunits->r.outscale << " " ;
  //Outflow << vel[3 *  part_ID + 0]*simunits->v.outscale << " " << vel[3 *  part_ID + 1]*simunits->v.outscale << " " << vel[3 *  part_ID + 2]*simunits->v.outscale << " " ;
  //Outflow << e_L[0] << " " << e_L[1] << " " << e_L[2] << " " <<  f_vel_dist << endl;
  //Outflow << pos[3 *  (part_ID+Npart) + 0]*simunits->r.outscale << " " << pos[3 *  (part_ID+Npart) + 1]*simunits->r.outscale << " " << pos[3 *  (part_ID+Npart) + 2]*simunits->r.outscale << " " ;
  //Outflow << vel[3 *  (part_ID+Npart) + 0]*simunits->v.outscale << " " << vel[3 *  (part_ID+Npart) + 1]*simunits->v.outscale << " " << vel[3 *  (part_ID+Npart) + 2]*simunits->v.outscale << " " ;
  //Outflow << e_L[0] << " " << e_L[1] << " " << e_L[2] << " " <<  f_vel_dist << endl;
  //Outflow.close();



  for (j=0; j<ndim; j++){
		vel[3 *  part_ID + j]	 		    += V_rot[j];
		vel[3 * (part_ID+Npart) + j] 	+= V_rot_sym[j];
	}

	return;
}


//=================================================================================================
//  Outflows::QuaternionMirror()
/// Rotation-Matrix to transform particle to angular-momentum-axis
//=================================================================================================
template <int ndim>
FLOAT* Outflow<ndim>::QuaternionMirror(FLOAT *axis, FLOAT *vec0)
{

	FLOAT		v[3]; 			// RoationAxis
	FLOAT		u[4];			  // Quaternion
	FLOAT		u1[4];			// inverse Quaternion
	FLOAT		vec1[4];		// intermediate result
	FLOAT		vec2[4];		// intermediate result
	FLOAT		vec3[4];		// intermediate result
	FLOAT		v_norm;			// norm of the rotation axis
	FLOAT		axis_norm;  // norm of the desired axis

	axis_norm 	= norm(axis);
	for (int k=0; k<3; k++) axis[k] /= axis_norm;

// Calculate the mirror axis
//-------------------------------------------------------------------------------------------
	v[0]	=	0.5*axis[0];
	v[1]	=	0.5*axis[1];
	v[2]	=	0.5*(axis[2]) + 0.5;

	v_norm 	= norm(v);
	for (int k=0; k<3; k++) v[k] /= v_norm;

// compute Quaternions u and u1
//-------------------------------------------------------------------------------------------
	vec1[0]	=	0.0;
	u[0]	=	0.0;
	u1[0]	=	0.0;
	for (int k=0; k<3; k++) u[k+1]		=	v[k];
	for (int k=0; k<3; k++) u1[k+1] 	= -u[k+1];
	for (int k=0; k<3; k++) vec1[k+1]	=	vec0[k];

// compute Hamilton product u*vec*u^-1
//-------------------------------------------------------------------------------------------
	Hprod(u, vec1, vec2);
	Hprod(vec2, u1, vec3);


	for (int k=0; k<3; k++)	vec0[k]  	=   vec3[k+1];
	return vec0;
}

//=================================================================================================
//  Outflows::Set_part()
/// Calls hydro->Create_New_Particle function
//=================================================================================================
template <int ndim>
void Outflow<ndim>::Set_part(
	int sink_ID,
	int Npart,
	int n,
	int level_step,
	int level_max,
	FLOAT v0,
	FLOAT current_time,
	FLOAT *pos,
	FLOAT *vel,
	FLOAT *pos_rev,
	FLOAT *vel_rev,
	Hydrodynamics<ndim> *hydro)
{

	int        i,j;
	FLOAT      r1[ndim];
	FLOAT      v1[ndim];
	FLOAT      r2[ndim];
	FLOAT      v2[ndim];
	FLOAT      v_abs2;
	FLOAT      u;
	Particle<ndim> *partPtr;

	for (j=0; j<2*Npart; j++){
		v_abs2 = 0.0;
		for (i=0; i<ndim; i++){
			r1[i] = pos[ndim*j+i];
			r2[i] = pos_rev[ndim*j+i];
			v1[i] = vel[ndim*j+i];
			v2[i] = vel_rev[ndim*j+i];
			v_abs2 += pow(vel[ndim*j+i],2);
			//assert(!(isnan(r1[i])));
			//assert(!(isnan(v1[i])));
			//assert(!(isnan(r2[i])));
			//assert(!(isnan(v2[i])));
		}

		u = 0.005*v_abs2;
		u = min(u, 1.0e8/simunits->u.outscale);



		partPtr = &hydro->CreateNewParticle(gas_type, dm, u, r1, v1, sim);
		partPtr->h = 6*sinks->sink[sink_ID].radius;
		partPtr->outflowFlag 	= current_time;
		cout << "outflow Time: " << partPtr->outflowFlag << endl;
		//cout << "Vsc1: " << v1[0]*simunits->v.outscale << " " << v1[1]*simunits->v.outscale << " " << v1[2]*simunits->v.outscale << " " << endl;
		//cout << " uscl: " <<  partPtr->u*simunits->u.outscale << "  u:  " << partPtr->u << " temp "
		// << " Et/Ek " << partPtr->u/(0.5*v_abs2) << endl;

		partPtr = &hydro->CreateNewParticle(gas_type, dm, u, r2, v2, sim);
		partPtr->h = 6*sinks->sink[sink_ID].radius;
		partPtr->outflowFlag 	= current_time;
		//cout << "Vsc2: " << v2[0]*simunits->v.outscale << " " << v2[1]*simunits->v.outscale << " " << v2[2]*simunits->v.outscale << " " << endl;
		//cout << " uscl: " <<  partPtr->u*simunits->u.outscale << "  u:  " << partPtr->u << " temp "
		// << " Et/Ek " << partPtr->u/(0.5*v_abs2) << endl;

	}
	return;
}

//=================================================================================================
//  Outflows::ComputeAngmom
/// Derive the break-up angular momentum of the star
/// and the amount each ejected particle should carray away.
//=================================================================================================
template <int ndim>
FLOAT Outflow<ndim>::ComputeAngmom(int SinkID, int Npart)
{
  FLOAT L_max;
  FLOAT D_L_max;
  FLOAT L_star;

  if(stellarEv !="none"){
    L_max     = sinks->sink[SinkID].M_star*sqrt(G_Const*sinks->sink[SinkID].M_star*sinks->sink[SinkID].stellarRadius);
  }else{
    L_max     = sinks->sink[SinkID].M_star*sqrt(G_Const*sinks->sink[SinkID].M_star*10*r_au / (r_pc*simunits->r.outscale));
  }
  L_star    = sqrt(dot_product(sinks->sink[SinkID].L_star, sinks->sink[SinkID].L_star));
  D_L_max   = max(L_star - f_angmom*L_max, 0.0);
  
  cout << "L_star " << L_star << " L_IAD " << sqrt(pow(sinks->sink[SinkID].L_IAD[0],2) + pow(sinks->sink[SinkID].L_IAD[1],2) + pow(sinks->sink[SinkID].L_IAD[2],2)) << " L_max  " << f_angmom*L_max << " " << max(L_star - f_angmom*L_max, 0.0) << endl;

  cout << "M_star " << sinks->sink[SinkID].M_star << " M_IAD " << sinks->sink[SinkID].star->m - sinks->sink[SinkID].M_star << endl;
  return 0.25*D_L_max;
}

//=================================================================================================
//  Outflows::Update_sink()
/// Subtract mass and angular momentum from the sink particle
//=================================================================================================
template <int ndim>
void Outflow<ndim>::Update_sink(int SinkID, int Npart, FLOAT *L_total)
{
	sinks->sink[SinkID].M_star -= 4.0 * Npart * dm;
	sinks->sink[SinkID].star->m -= 4.0 * Npart * dm;
	for (int k=0; k<ndim; k++) sinks->sink[SinkID].angmom[k] -= L_total[k];
	for (int k=0; k<ndim; k++) sinks->sink[SinkID].L_star[k] -= L_total[k];
	for (int k=0; k<ndim; k++) sinks->sink[SinkID].L_IAD[k] = sinks->sink[SinkID].angmom[k] - sinks->sink[SinkID].L_star[k];

	return;
}

//=================================================================================================
//  Outflows::CallOutput()
/// Prints inportant information about the outflow
//=================================================================================================
template <int ndim>
void Outflow<ndim>::CallOutput(int SinkID, int Npart, FLOAT *L_total, FLOAT *e_L, FLOAT v0)
{
	cout << "v0 		: " << v0*simunits->v.outscale													<< endl;
	cout << "L_sink 	: " << sinks->sink[SinkID].angmom[0]*simunits->angmom.outscale 	<< " "
							<< sinks->sink[SinkID].angmom[1]*simunits->angmom.outscale 	<< " "
							<< sinks->sink[SinkID].angmom[2]*simunits->angmom.outscale 	<< endl;
	cout << "L_IAD  	: " << sinks->sink[SinkID].L_IAD[0]*simunits->angmom.outscale  	<< " "
							<< sinks->sink[SinkID].L_IAD[1]*simunits->angmom.outscale  	<< " "
							<< sinks->sink[SinkID].L_IAD[2]*simunits->angmom.outscale  	<< endl;
	cout << "L_star 	: " << sinks->sink[SinkID].L_star[0]*simunits->angmom.outscale 	<< " "
							<< sinks->sink[SinkID].L_star[1]*simunits->angmom.outscale 	<< " "
							<< sinks->sink[SinkID].L_star[2]*simunits->angmom.outscale 	<< endl;
	cout << "e_L 		: " <<  e_L[0] << " " << e_L[1] << " " << e_L[2] 				<< endl;
	cout << "L_out		: "	<< L_total[0]*simunits->angmom.outscale						<< " "
							<< L_total[1]*simunits->angmom.outscale						<< " "
							<< L_total[2]*simunits->angmom.outscale						<< endl;
	cout << "|L_out|	: "	<< sqrt(pow(L_total[0],2)
								+pow(L_total[1],2)
								+pow(L_total[2],2))*simunits->angmom.outscale			<< endl;
	cout << "|L_IAD|	: "	<< sqrt(pow(sinks->sink[SinkID].L_IAD[0],2)
								+pow(sinks->sink[SinkID].L_IAD[1],2)
								+pow(sinks->sink[SinkID].L_IAD[2],2))
								*simunits->angmom.outscale 								<< endl;
	cout << "|L_star|	: "	<< sqrt(pow(sinks->sink[SinkID].L_star[0],2)
								+pow(sinks->sink[SinkID].L_star[1],2)
								+pow(sinks->sink[SinkID].L_star[2],2))
								*simunits->angmom.outscale 								<< endl;
	cout << "d_L_tot   	: " << sinks->sink[SinkID].dL_Outflow * 4*Npart
								*simunits->angmom.outscale								<< endl;
	return;
}

//=================================================================================================
//  Outflows::dot_product()
/// Simple scalar-product between two vectors
//=================================================================================================
template <int ndim>
FLOAT Outflow<ndim>::dot_product(FLOAT *vec1, FLOAT *vec2)
{
	int 		j;
	FLOAT 		product;

	product = 0.0;
	for (j=0; j<ndim; j++){
		product += vec1[j] * vec2[j];
	}
	return product;
}

//=================================================================================================
//  Outflows::cross_product()
/// Simple cross-product between two vectors
//=================================================================================================
template <int ndim>
FLOAT* Outflow<ndim>::cross_product(FLOAT *vec1, FLOAT *vec2, FLOAT *vec3)
{
	vec3[0] = (vec1[1]*vec2[2] - vec1[2]*vec2[1]);
	vec3[1] = (vec1[2]*vec2[0] - vec1[0]*vec2[2]);
	vec3[2] = (vec1[0]*vec2[1] - vec1[1]*vec2[0]);

	return vec3;
}

//=================================================================================================
//  Outflows::Hprod()
/// Hamilton Product of Quaternions
//=================================================================================================
template <int ndim>
FLOAT* Outflow<ndim>::Hprod(FLOAT *vec1, FLOAT *vec2, FLOAT *vec3)
{

	vec3[0] = vec1[0]*vec2[0] - vec1[1]*vec2[1] - vec1[2]*vec2[2] - vec1[3]*vec2[3];
	vec3[1] = vec1[0]*vec2[1] + vec1[1]*vec2[0] + vec1[2]*vec2[3] - vec1[3]*vec2[2];
	vec3[2] = vec1[0]*vec2[2] - vec1[1]*vec2[3] + vec1[2]*vec2[0] + vec1[3]*vec2[1];
	vec3[3] = vec1[0]*vec2[3] + vec1[1]*vec2[2] - vec1[2]*vec2[1] + vec1[3]*vec2[0];

	return vec3;
}

//=================================================================================================
//  Outflows::norm()
/// Calculates the absolute value of a vector
//=================================================================================================
template <int ndim>
FLOAT Outflow<ndim>::norm(FLOAT *vec)
{
	int 		j;
	FLOAT 		norm;

	norm = 0.0;
	for (j=0; j<ndim; j++){
		norm += vec[j] * vec[j];
	}
	return sqrt(norm);
}



template class Outflow<1>;
template class Outflow<2>;
template class Outflow<3>;









