//=================================================================================================
//  SphSimulation.cpp
//  Contains all main functions controlling SPH simulation work-flow.
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
#include <iomanip>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include "Precision.h"
#include "CodeTiming.h"
#include "Exception.h"
#include "Debug.h"
#include "Ic.h"
#include "InlineFuncs.h"
#include "Simulation.h"
#include "Parameters.h"
#include "Nbody.h"
#include "RandomNumber.h"
#include "Sph.h"
#include "RiemannSolver.h"
#include "Ghosts.h"
#include "Sinks.h"
#include "StellarEvolution.h"
using namespace std;



// Create template class instances of the main SphSimulation object for
// each dimension used (1, 2 and 3)
template class SphSimulation<1>;
template class SphSimulation<2>;
template class SphSimulation<3>;



//=================================================================================================
//  SphSimulation::ProcessParameters
/// Process all the options chosen in the parameters file, setting various
/// simulation variables and creating important simulation objects.
/// SPH specific version
//=================================================================================================
template <int ndim>
void SphSimulation<ndim>::ProcessParameters(void)
{
  // Local references to parameter variables for brevity
  map<string, int> &intparams = simparams->intparams;
  map<string, double> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;
  string sim = stringparams["sim"];
  string gas_eos = stringparams["gas_eos"];
  string gas_radiation = stringparams["radiation"];

  debug2("[SphSimulation::ProcessParameters]");

  // Common set-up
  Simulation<ndim>::ProcessParameters();


  // Set-up main SPH objects depending on which SPH algorithm we are using
  ProcessSphParameters();
  hydro = sph;
  neib = sphneib ;

  // Process all N-body parameters and set-up main N-body objects
  this->ProcessNbodyParameters();


  // Set pointers to external potential field object
  sph->extpot = extpot;
  nbody->extpot = extpot;

  // Set eos pointer to nbody
  sph->eos->set_nbody_data(nbody);

  // Set all other SPH parameter variables
  sph->Nhydromax       = intparams["Nhydromax"];
  sph->create_sinks    = intparams["create_sinks"];
  sph->fixed_sink_mass = intparams["fixed_sink_mass"];
  sph->msink_fixed     = floatparams["m1"];
  sph->alpha_visc_min  = floatparams["alpha_visc_min"];


  // Set important variables for N-body objects
  nbody->Nstarmax       = intparams["Nstarmax"];
  nbody_single_timestep = intparams["nbody_single_timestep"];
  nbodytree.gpehard     = floatparams["gpehard"];
  nbodytree.gpesoft     = floatparams["gpesoft"];
  //nbody->perturbers     = intparams["perturbers"];
  //if (intparams["sub_systems"] == 1) subsystem->perturbers = intparams["perturbers"];


  // Sink particles
  //-----------------------------------------------------------------------------------------------
  sinks = new Sinks<ndim>(sphneib);
  sink_particles             = intparams["sink_particles"];
  sinks->sink_particles      = intparams["sink_particles"];
  sinks->create_sinks        = intparams["create_sinks"];
  sinks->smooth_accretion    = intparams["smooth_accretion"];
  sinks->Nsinkfixed          = intparams["Nsinkfixed"];
  sinks->alpha_ss            = floatparams["alpha_ss"];
  sinks->smooth_accrete_frac = floatparams["smooth_accrete_frac"];
  sinks->smooth_accrete_dt   = floatparams["smooth_accrete_dt"];
  sinks->sink_radius_mode    = stringparams["sink_radius_mode"];
  sinks->rho_sink            = floatparams["rho_sink"];
  sinks->rho_sink            /= (simunits.rho.outscale*simunits.rho.outcgs);

  if (sinks->sink_radius_mode == "fixed") {
    sinks->sink_radius = floatparams["sink_radius"]/simunits.r.outscale;
  }
  else {
    sinks->sink_radius = floatparams["sink_radius"];
  }

  // Sanity-check for various sink particle values
  if (intparams["sink_particles"] == 1 &&
      (stringparams["nbody"] != "lfkdk" && stringparams["nbody"] != "lfdkd")) {
    string message = "Invalid parameter : nbody must use lfkdk or lfdkd when "
      "using accreting sink particles";
    ExceptionHandler::getIstance().raise(message);
  }

#if defined MPI_PARALLEL
  sinks->SetMpiControl(mpicontrol);
  if (stringparams["out_file_form"]=="sf") {
    string message = "The sf format is not supported with MPI! Use the column "
        "or (better) the su format";
    ExceptionHandler::getIstance().raise(message);
  }
#endif

  // Supernova feedback
  //-----------------------------------------------------------------------------------------------
  if (stringparams["supernova_feedback"] == "none") {
    snDriver = new NullSupernovaDriver<ndim>(this);
  }
  else if (stringparams["supernova_feedback"] == "single") {
    snDriver = new SedovTestDriver<ndim>(this, simparams, simunits);
  }
  else if (stringparams["supernova_feedback"] == "random") {
    snDriver = new RandomSedovTestDriver<ndim>(this, simparams, simunits, simbox);
  }
  else if (stringparams["supernova_feedback"] == "silcc") {
    snDriver = new SilccSupernovaDriver<ndim>(this, simparams, simunits);
  }
  else {
    string message = "Unrecognised parameter :  = supernova_feedback"
      + simparams->stringparams["supernova_feedback"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Accretion Feedback
  //-----------------------------------------------------------------------------------------------
	if (stringparams["acc_fb"] != "none") {
		if (stringparams["acc_fb"] == "episodic") {
			accfb = new EpisodicAFB<ndim>(&simunits, this, stringparams["stellarEv"], floatparams["AccFBminMass"], floatparams["dmdt_BGR"],
				floatparams["alpha_MRI"], floatparams["gamma_eos"], floatparams["mu_bar"],
				floatparams["temp0"], floatparams["f_angmom"]);
			if (intparams["outflow_fb"] == 1) {
				outf = new Outflow<ndim>(&simunits, this, randnumb, simparams->stringparams["stellarEv"], intparams["Nout"],
          floatparams["R_launch"], floatparams["f_vel"], floatparams["op_angle"],
          floatparams["theta0"], floatparams["f_angmom"], floatparams["f_eject"], floatparams["r_eject"]);
			}
		}
  	else if (stringparams["acc_fb"] == "cont") {
    	accfb = new ContinuousAFB<ndim>(&simunits, this, stringparams["stellarEv"], floatparams["gamma_eos"], floatparams["mu_bar"],
				floatparams["temp0"], floatparams["dmdt_max"], floatparams["f_angmom"], floatparams["AccFBminMass"]);
			if (intparams["outflow_fb"] == 1) {
				outf = new Outflow<ndim>(&simunits, this, randnumb, simparams->stringparams["stellarEv"], intparams["Nout"],
        floatparams["R_launch"], floatparams["f_vel"], floatparams["op_angle"],
        floatparams["theta0"], floatparams["f_angmom"], floatparams["f_eject"], floatparams["r_eject"]);
			}
		}
    else{
      cout << "Attetion : No/wrong accretion feedback mode " << endl;
    }
	}

  if (simparams->stringparams["stellarEv"] != "none"){
  	stellarEv = new StellarEvolution<ndim>(&simunits, this, floatparams["AccFBminMass"], stringparams["acc_fb"] );
  }

  // Set pointers to timing object
  nbody->timing   = timing;
  if (sim == "sph" || sim == "gradhsph" || sim == "sm2012sph") {
    sinks->timing     = timing;
    hydro->timing     = timing;
    hydroint->timing  = timing;
    uint->timing      = timing;
    radiation->timing = timing;
    sphneib->SetTimingObject(timing);
  }


  // Flag that we've processed all parameters already
  ParametersProcessed = true;


  return;
}



//=================================================================================================
//  SphSimulation::PostInitialConditionsSetup
/// Call routines for calculating all initial SPH and N-body quantities
/// once initial conditions have been set-up.
//=================================================================================================
template <int ndim>
void SphSimulation<ndim>::PostInitialConditionsSetup(void)
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  FLOAT adot[ndim];                    // Dummy adot array

  debug2("[SphSimulation::PostInitialConditionsSetup]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_POST_IC_SETUP");

  sph->DeleteDeadParticles();

  // Set iorig
  if (rank == 0) {
    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).iorig = i;
  }

  // Copy information from the stars to the sinks
  if (simparams->intparams["sink_particles"]==1) {
    sinks->Nsink = nbody->Nstar;
    sinks->AllocateMemory(sinks->Nsink);
    for (int i=0; i<sinks->Nsink; i++) {
      sinks->sink[i].star   = &(nbody->stardata[i]);
      sinks->sink[i].istar  = i;
      sinks->sink[i].radius = hydro->kernp->kernrange*nbody->stardata[i].h;
      //sinks->sink[i].radius = simparams->floatparams["sink_radius"];
    }
  }

	uint->sinks = sinks;

  // Perform initial MPI decomposition
  //-----------------------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
  mpicontrol->CreateInitialDomainDecomposition(sph, nbody, simparams,
                                               this->initial_h_provided);
  this->AllocateParticleMemory();
#endif


  // Set time variables here (for now)
  nresync = 0;
  n = 0;

  // Set initial smoothing lengths and create initial ghost particles
  //-----------------------------------------------------------------------------------------------
  sph->Nghost = 0;
  sph->Ntot = sph->Nhydro;
  for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).flags.set(active);

  // Set initial artificial viscosity alpha values
  if (simparams->stringparams["time_dependent_avisc"] == "none") {
    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).alpha = sph->alpha_visc;
  }
  else {
    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).alpha = sph->alpha_visc_min;
  }

  // Compute instantaneous mean mass (used for smooth sink accretion)
  sph->mmean = (FLOAT) 0.0;
  for (i=0; i<sph->Nhydro; i++) sph->mmean += sph->GetSphParticlePointer(i).m;
  sph->mmean /= (FLOAT) sph->Nhydro;

  // Compute minimum smoothing length of sinks
  sph->hmin_sink = big_number;
  for (i=0; i<sinks->Nsink; i++) {
    sph->hmin_sink = min(sph->hmin_sink, (FLOAT) sinks->sink[i].star->h);
  }

  // If the smoothing lengths have not been provided beforehand, then
  // calculate the initial values here
  //sphneib->neibcheck = false;
  if (!this->initial_h_provided) {
    sph->InitialSmoothingLengthGuess();
    sphneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, timestep, sph);
#ifdef MPI_PARALLEL
  sphneib->InitialiseCellWorkCounters();
#endif
    sphneib->UpdateAllSphProperties(sph, nbody, simbox);
  }
  else {
    sphneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, timestep, sph);
#ifdef MPI_PARALLEL
  sphneib->InitialiseCellWorkCounters();
#endif
  }


#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro, sph, sph->kernp,false, false);
#endif

  // Search ghost particles
  sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, sph);
#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro+sph->NPeriodicGhost, sph, sph->kernp,true,false);
  for (int i=0; i<sph->Nhydro+sph->NPeriodicGhost; i++) {
    SphParticle<ndim>& parti = sph->GetSphParticlePointer(i);
    parti.hrangesqd = sph->kernp->kernrangesqd*parti.h*parti.h;
  }
  MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep,  sph);
#endif

  // Zero accelerations
  for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).flags.set(active);

  // Update neighbour tree
  rebuild_tree = true;
  sphneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, timestep, sph);
#ifdef MPI_PARALLEL
  sphneib->InitialiseCellWorkCounters();
#endif
  level_step = 1;


  // For Eigenvalue MAC, need non-zero values
  for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).gpot = big_number;

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph, nbody, simbox);

#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro, sph, sph->kernp,false,true);
#endif

  // Update neighbour tree
  rebuild_tree = true;
  sphneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, timestep, sph);
#ifdef MPI_PARALLEL
  sphneib->InitialiseCellWorkCounters();
#endif
  sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, sph);

#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro + sph->NPeriodicGhost, sph, sph->kernp,true,true);
  MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, sph);
#endif

    // Communicate pruned trees for MPI
#ifdef MPI_PARALLEL
  sphneib->BuildPrunedTree(rank, simbox, mpicontrol->mpinode, sph);
#endif


  // Compute all initial N-body terms
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<nbody->Nstar; i++) {
    for (k=0; k<ndim; k++) nbody->stardata[i].a[k]     = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) nbody->stardata[i].adot[k]  = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) nbody->stardata[i].a2dot[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) nbody->stardata[i].a3dot[k] = (FLOAT) 0.0;
    nbody->stardata[i].gpot   = (FLOAT) 0.0;
    nbody->stardata[i].gpe    = (FLOAT) 0.0;
    nbody->stardata[i].tlast  = t;
    nbody->stardata[i].flags.set(active);
    nbody->stardata[i].level  = 0;
    nbody->stardata[i].nstep  = 0;
    nbody->stardata[i].nlast  = 0;
    nbody->nbodydata[i]       = &(nbody->stardata[i]);
  }
  nbody->Nnbody = nbody->Nstar;


  // Read-in stellar properties table here
  nbody->LoadStellarPropertiesTable(&simunits);
  nbody->UpdateStellarProperties();


  // Compute all initial SPH force terms
  // We will need to iterate if we are going to use a relative opening criterion
  //-----------------------------------------------------------------------------------------------
  MAC_Type mac_type = sphneib->GetOpeningCriterion() ;
  int n_iter = 1 + (sph->self_gravity == 1 && mac_type != geometric) ;
  for (int iter=0; iter < n_iter; iter++) {

    if (iter==0 && mac_type != geometric)
      sphneib->SetOpeningCriterion(geometric);
    else
      sphneib->SetOpeningCriterion(mac_type) ;

    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).flags.set(active);


    // Copy all other data from real SPH particles to ghosts
    LocalGhosts->CopyHydroDataToGhosts(simbox, sph);

    sphneib->BuildTree(true, 0, ntreebuildstep, ntreestockstep, timestep, sph);
#ifdef MPI_PARALLEL
    sphneib->InitialiseCellWorkCounters();
#endif
    sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, sph);
    sphneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, sph);
#ifdef MPI_PARALLEL
    mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro + sph->NPeriodicGhost, sph, sph->kernp,true,true);
    MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, sph);
    sphneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, sph);
#endif

    // Zero accelerations (here for now)
    for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      part.level     = 0;
      part.levelneib = 0;
      part.nstep     = 0;
      part.nlast     = 0;
      part.tdav			 = (FLOAT) 0.0;
      part.div_v     = (FLOAT) 0.0;
      part.dudt      = (FLOAT) 0.0;
      part.gpot      = (FLOAT) 0.0;
      part.mu_bar    = (FLOAT) simparams->floatparams["mu_bar"];
      for (k=0; k<ndim; k++) part.a[k] = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) part.atree[k] = (FLOAT) 0.0;
    }

    // Calculate all SPH properties
    if (iter == 0)
      sphneib->UpdateAllSphProperties(sph, nbody, simbox);


#ifdef MPI_PARALLEL
    if (sph->self_gravity == 1) {
      sphneib->UpdateGravityExportList(rank, sph, nbody, simbox, ewald);
    }
    else {
      sphneib->UpdateHydroExportList(rank, sph, nbody, simbox);
    }

    mpicontrol->ExportParticlesBeforeForceLoop(sph);
#endif

    if (iter == 0) {
      // Update the radiation field
      for (int jj=0; jj<10; jj++) {
        radiation->UpdateRadiationField(sph->Nhydro, nbody->Nnbody, sinks->Nsink,
            sph->GetSphParticleArray(), nbody->nbodydata, sinks->sink);
      }

      // Update thermal properties (if radiation field has altered them)
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        sph->ComputeThermalProperties(part);
      }
    }

    // Calculate SPH gravity and hydro forces, depending on which are activated
    if (sph->self_gravity == 1) {
      sphneib->UpdateAllSphForces(sph, nbody, simbox, ewald);
    }
    else if (sph->hydro_forces == 1) {
      sphneib->UpdateAllSphHydroForces(sph, nbody, simbox);
    }
    else{
      ExceptionHandler::getIstance().raise("Error: No forces included in simulation");
    }

#if defined MPI_PARALLEL
    mpicontrol->GetExportedParticlesAccelerations(sph);
#endif
  } // End of iteration loop.

  // Add external potential for all active SPH particles
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    sph->extpot->AddExternalPotential(part.r, part.v, part.a, adot, part.gpot);
  }

  // Set initial accelerations
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    for (k=0; k<ndim; k++) part.r0[k] = part.r[k];
    for (k=0; k<ndim; k++) part.v0[k] = part.v[k];
    for (k=0; k<ndim; k++) part.a0[k] = part.a[k];
    part.dudt0 = part.dudt;

    part.flags.unset(active);
  }


  // Compute initial N-body forces
  //-----------------------------------------------------------------------------------------------
  if (sph->self_gravity == 1 && sph->Nhydro > 0) {
    sphneib->UpdateAllStarGasForces(sph, nbody, simbox, ewald);
#if defined MPI_PARALLEL
    // We need to sum up the contributions from the different domains
    mpicontrol->ComputeTotalStarGasForces(nbody);
#endif
  }

  if (nbody->nbody_softening == 1) {
    nbody->CalculateDirectSmoothedGravForces(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
  }
  else {
    nbody->CalculateDirectGravForces(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
  }
  nbody->CalculateAllStartupQuantities(nbody->Nnbody, nbody->nbodydata, simbox, ewald);

  for (i=0; i<nbody->Nnbody; i++) {
    if (nbody->nbodydata[i]->flags.check(active)) {
      nbody->extpot->AddExternalPotential(nbody->nbodydata[i]->r, nbody->nbodydata[i]->v,
                                          nbody->nbodydata[i]->a, nbody->nbodydata[i]->adot,
                                          nbody->nbodydata[i]->gpot);
    }
  }

  // Compute the dust forces if present.
  if (sphdust != NULL){

    // Do a first guess of the drag forces
    for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      part.flags.set(active);
    }

    // Copy properties from original particles to ghost particles
    LocalGhosts->CopyHydroDataToGhosts(simbox, sph);
 #ifdef MPI_PARALLEL
    MpiGhosts->CopyHydroDataToGhosts(simbox, sph);
 #endif


    sphdust->UpdateAllDragForces(sph) ;
  }


  // Compute timesteps for all particles and use it to compute the time-averaged initial
  // drag force.
  if (Nlevels == 1) this->ComputeGlobalTimestep();
  else this->ComputeBlockTimesteps();

  // Compute the time-averaged initial dust forces.
  if (sphdust != NULL) {

    // Set the acceleration to the pre-drag value
    for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      for (k=0; k <ndim; k++) part.a[k] = part.a0[k] ;
      part.dudt = part.dudt0;
    }

    sphdust->UpdateAllDragForces(sph) ;
  }

  // Set particle values for initial step (e.g. r0, v0, a0, u0)
  uint->EndTimestep(n, t, timestep, sph);
  hydroint->EndTimestep(n,t, timestep, sph);
  nbody->EndTimestep(n, nbody->Nstar, t, timestep, nbody->nbodydata);

  this->CalculateDiagnostics();
  this->diag0 = this->diag;
  this->setup = true;



  return;
}



//=================================================================================================
//  SphSimulation::MainLoop
/// Main SPH simulation integration loop.
//=================================================================================================
template <int ndim>
void SphSimulation<ndim>::MainLoop(void)
{
  int activecount = 0;                 // Flag if we need to recompute particles
  int i;                               // Particle loop counter
  int it;                              // Time-symmetric iteration counter
  int k;                               // Dimension counter
  FLOAT adot[ndim];                    // Dummy adot variable
  FLOAT tghost;                        // Approx. ghost particle lifetime

  debug2("[SphSimulation::MainLoop]");

  // Advance time variables
  n = n + 1;
  Nsteps = Nsteps + 1;
  t = t + timestep;
  if (n == nresync) Nblocksteps = Nblocksteps + 1;
  if (n%integration_step == 0) Nfullsteps = Nfullsteps + 1;

  // Advance SPH and N-body particles' positions and velocities
  uint->EnergyIntegration(n, (FLOAT) t, (FLOAT) timestep, sph);
  hydroint->AdvanceParticles(n, (FLOAT) t, (FLOAT) timestep, sph);
  nbody->AdvanceParticles(n, nbody->Nnbody, t, timestep, nbody->nbodydata);


  // Add any new particles into the simulation here (e.g. Supernova, wind feedback, etc..).
  //-----------------------------------------------------------------------------------------------

	//cout << "n " << n << " pow " << pow(2,level_step - level_max) << " level_step " << level_step << " level_max " << level_max << " nresync " << nresync << " integration_step " << integration_step << " Nlevels " << Nlevels <<  endl;
   //if (n%(int) pow(2,level_step - level_max) == 0) {

	if (n%(int) pow(2,level_step - level_max) == 0) {
  	if (simparams->stringparams["supernova_feedback"] != "none") {
    	snDriver->Update(n, level_step, level_max, t, hydro, sphneib, randnumb);
  	}
	}
	if (n == nresync) outflowTiming = false;
  if (n == nresync || outflowTiming) {
  	if (simparams->stringparams["acc_fb"] != "none" and sinks->Nsink > 0) {
			accfb->Luminosity(hydro, sinks, simparams->floatparams["f_acc"]);
  	}
		if (simparams->intparams["outflow_fb"] == 1 and sinks->Nsink > 0) {
			outf->Check_Outflow(sinks, hydro, t, n, level_step, level_max);
		}
    if (sph->newParticles){
			outflowTiming = true;
    	rebuild_tree = true;
    	sphneib->UpdateAllSphProperties(sph, nbody, simbox);
  		//if (Nlevels == 1) this->ComputeGlobalTimestep();
  		//else this->ComputeBlockTimesteps();
  		sph->newParticles = false;
  		cout << "new SPH part" << endl;
  	}
  }

  // Check all boundary conditions
  // (DAVID : create an analagous of this function for N-body)
  hydroint->CheckBoundaries(simbox,sph);


  // Perform the load-balancing step for MPI simulations.  Need to stock local and pruned trees,
  // then compute the new load-balanced MPI domains and finally transfer the
  // particles to the new domains.
  //-----------------------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
  if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
	// Horrible hack in order NOT to trigger a full tree rebuild
	int Nstepsaux=Nsteps;
	if (Nstepsaux%2==0) Nstepsaux++;
	sphneib->BuildTree(rebuild_tree,Nstepsaux,2, ntreestockstep,timestep,sph);
	if (rebuild_tree) {
		  sphneib->BuildPrunedTree(rank, simbox, mpicontrol->mpinode, sph);
	}
	else {
		sphneib->StockPrunedTree(rank, sph);
	}
    mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro, sph, sph->kernp,false,true);
    mpicontrol->LoadBalancing(sph, nbody);
  }
#endif


  // Rebuild or update local neighbour and gravity tree
  sphneib->BuildTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep, timestep, sph);

#ifdef MPI_PARALLEL
  sphneib->InitialiseCellWorkCounters();
#endif

  // Search for new ghost particles and create on local processor
  //-----------------------------------------------------------------------------------------------
  //tghost = timestep*(FLOAT)(ntreebuildstep - 1);
  tghost = 0;
  sphneib->SearchBoundaryGhostParticles(tghost, simbox, sph);
  sphneib->BuildGhostTree(true, Nsteps, ntreebuildstep, ntreestockstep, timestep, sph);


  // Re-build and communicate the new pruned trees (since the trees will necessarily change
  // once there has been communication of particles to new domains)
#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro + sph->NPeriodicGhost, sph, sph->kernp,true,true);
{
CodeTiming::BlockTimer timer = timing->StartNewTimer("SEARCH_GHOSTS_MPI");
  MpiGhosts->SearchGhostParticles(tghost, simbox, sph);
}
  sphneib->BuildMpiGhostTree(true, Nsteps, ntreebuildstep, ntreestockstep, timestep, sph);
#endif


  // Iterate if we need to immediately change SPH particle timesteps
  // (e.g. due to feedback, or sudden change in neighbour timesteps)
  //-----------------------------------------------------------------------------------------------
  do {

    // Update cells containing active particles
    if (activecount > 0) sphneib->UpdateActiveParticleCounters(sph);

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph, nbody, simbox);

    // Zero accelerations (here for now)
    sph->ZeroAccelerations();

    // Update the radiation field
    if (simparams->stringparams["radiation"] != "none" && (Nsteps%nradstep == 0 || recomputeRadiation)) {
      radiation->UpdateRadiationField(sph->Nhydro, nbody->Nnbody, sinks->Nsink,
                                      sph->GetSphParticleArray(), nbody->nbodydata, sinks->sink);
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        sph->ComputeThermalProperties(part);
      }
    }

    if (simparams->intparams["acc_fb"] == 1 and sinks->Nsink > 0) {
    		accfb->updateAccFBTemp(sph->Nhydro, hydro, sinks);
    }

    // Calculate gravitational forces from other distant MPI nodes.
    // Also determines particles that must be exported to other nodes
    // if too close to the domain boundaries
#ifdef MPI_PARALLEL
    // Pruned trees are used only to compute which particles to export
    // Therefore we don't need to update them at the start of the loop, and we can do it soon before we need them
    if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
      sphneib->BuildPrunedTree(rank, simbox, mpicontrol->mpinode, sph);
    }
    else {
      sphneib->StockPrunedTree(rank, sph);
    }

    if (sph->self_gravity == 1) {
      sphneib->UpdateGravityExportList(rank, sph, nbody, simbox, ewald);
    }
    else {
      sphneib->UpdateHydroExportList(rank, sph, nbody, simbox);
    }

    // If active particles need forces from other domains, export particles
    mpicontrol->ExportParticlesBeforeForceLoop(sph);
#endif
    // Calculate SPH gravity and hydro forces, depending on which are activated
    if (sph->self_gravity == 1) {
      sphneib->UpdateAllSphForces(sph, nbody, simbox, ewald);
    }
    else if (sph->hydro_forces == 1) {
      sphneib->UpdateAllSphHydroForces(sph, nbody, simbox);
    }


    // Add external potential for all active SPH particles
    if (simparams->stringparams["external_potential"] != "none") {
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.check(active)) {
          sph->extpot->AddExternalPotential(part.r, part.v, part.a, adot, part.gpot);
        }
      }
    }


    // Checking if acceleration or other values are invalid
#ifndef NDEBUG
    for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      if (part.flags.check(active)) {
        for (k=0; k<ndim; k++) assert(isfinite(part.r[k]));
        for (k=0; k<ndim; k++) assert(isfinite(part.v[k]));
        for (k=0; k<ndim; k++) assert(isfinite(part.a[k]));
        assert(isfinite(part.gpot));
      }
    }
#endif

#if defined MPI_PARALLEL
    mpicontrol->GetExportedParticlesAccelerations(sph);
#endif

    // Check if all neighbouring timesteps are acceptable.  If not, then set any
    // invalid particles to active to recompute forces immediately.
    if (Nlevels > 1) {
      activecount = hydroint->CheckTimesteps(level_diff_max, level_step, n, timestep, sph);
    }
    else {
      activecount = 0;
    }

#if defined MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &activecount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif


  } while (activecount > 0);
  //-----------------------------------------------------------------------------------------------


  // Check that we have sensible smoothing lengths
  if (periodicBoundaries) {
    double hmax = sphneib->GetMaximumSmoothingLength() ;
    hmax *= sph->kernp->kernrange ;
    for (int k=0; k<ndim; k++) {
      if (simbox.half[k] < 2*hmax) {
        string message = "Error: Smoothing length too large, self-interaction will occur";
        ExceptionHandler::getIstance().raise(message);
      }
    }
  }

  // Iterate for P(EC)^n schemes for N-body particles
  //-----------------------------------------------------------------------------------------------
  for (it=0; it<nbody->Npec; it++) {

    // Zero all acceleration terms
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->flags.check(active)) {
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a[k]     = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->adot[k]  = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a2dot[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a3dot[k] = (FLOAT) 0.0;
        nbody->nbodydata[i]->gpot = (FLOAT) 0.0;
        nbody->nbodydata[i]->gpe = (FLOAT) 0.0;
      }
    }
    if (sph->self_gravity == 1) {
      if (sph->Nhydro>0) sphneib->UpdateAllStarGasForces(sph, nbody, simbox, ewald);
#if defined MPI_PARALLEL
      // We need to sum up the contributions from the different domains
      mpicontrol->ComputeTotalStarGasForces(nbody);
#endif
    }

    // Calculate forces, force derivatives etc.., for active stars/systems
    if (nbody->nbody_softening == 1) {
      nbody->CalculateDirectSmoothedGravForces(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
    }
    else {
      nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata, simbox, ewald);
    }

    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->flags.check(active)) {
        nbody->extpot->AddExternalPotential(nbody->nbodydata[i]->r, nbody->nbodydata[i]->v,
                                            nbody->nbodydata[i]->a, nbody->nbodydata[i]->adot,
                                            nbody->nbodydata[i]->gpot);
      }
    }

    // Calculate correction step for all stars at end of step, except the
    // final iteration (since correction is computed in EndStep also).
    //if (it < nbody->Npec - 1)
    nbody->CorrectionTerms(n, nbody->Nnbody, t, timestep, nbody->nbodydata);

  }
  //-----------------------------------------------------------------------------------------------

  // Search for new sink particles (if activated) and accrete to existing sinks
  if (sink_particles == 1) {
    if (sinks->create_sinks == 1 && (rebuild_tree || Nfullsteps%ntreebuildstep == 0)) {
      sinks->SearchForNewSinkParticles(n, t ,sph, nbody);
    }
    if (sinks->Nsink > 0) {
      sph->mmean = (FLOAT) 0.0;
      for (i=0; i<sph->Nhydro; i++) sph->mmean += sph->GetSphParticlePointer(i).m;
      sph->mmean /= (FLOAT) sph->Nhydro;
      sph->hmin_sink = big_number;
      for (i=0; i<sinks->Nsink; i++) {
        sph->hmin_sink = min(sph->hmin_sink, (FLOAT) sinks->sink[i].star->h);
      }

      sinks->AccreteMassToSinks(n, timestep, sph, nbody);
      nbody->UpdateStellarProperties();

  		if (simparams->stringparams["stellarEv"] != "none"){
  			stellarEv->StellarEvolutionMain(sinks);
  		}
      if (extra_sink_output) WriteExtraSinkOutput();
    }
  }

  // Compute timesteps for all particles (for next step, needed for dust accelerations.)
  if (Nlevels == 1) this->ComputeGlobalTimestep();
  else this->ComputeBlockTimesteps();

  // Compute the dust forces if present.
  if (sphdust != NULL) {
    // Need to reset active particles:
    hydroint->SetActiveParticles(n, sph) ;
    sphneib->UpdateActiveParticleCounters(sph);

    // Copy properties from original particles to ghost particles
    LocalGhosts->CopyHydroDataToGhosts(simbox, sph);
#ifdef MPI_PARALLEL
    MpiGhosts->CopyHydroDataToGhosts(simbox, sph);
#endif

    sphdust->UpdateAllDragForces(sph) ;

    // Unset active particles
    for (i=0; i<sph->Nhydro; i++) {
      sph->GetSphParticlePointer(i).flags.unset(active);
    }
  }



  // Reset flags in preparation for next timestep
  rebuild_tree        = false;
  recomputeRadiation  = false;


  // End-step terms for all SPH particles
  if (sph->Nhydro > 0) {
    uint->EndTimestep(n, (FLOAT) t, (FLOAT) timestep, sph);
    hydroint->EndTimestep(n, (FLOAT) t, (FLOAT) timestep, sph);
  }

  // End-step terms for all star particles
  if (nbody->Nstar > 0) nbody->EndTimestep(n, nbody->Nnbody, t, timestep, nbody->nbodydata);


  return;
}





//=================================================================================================
//  SphSimulation::WriteExtraSinkOutput
/// For any simulations loaded into memory via a snapshot file, all particle
/// variables are converted into dimensionless code units here.
//=================================================================================================
template <int ndim>
void SphSimulation<ndim>::WriteExtraSinkOutput(void)
{
  int k;                               // ..
  int s;                               // ..
  string filename;                     // Output snapshot filename
  string nostring;                     // String of number of snapshots
  stringstream ss;                     // Stream object for preparing filename
  ofstream outfile;                    // Stream of restart file


  //-----------------------------------------------------------------------------------------------
  for (s=0; s<sinks->Nsink; s++) {

    SinkParticle<ndim> &sink = sinks->sink[s];
    nostring = "";
    ss << setfill('0') << setw(5) << s;
    nostring = ss.str();
    filename = run_id + ".sink." + nostring;
    ss.str(std::string());

    outfile.open(filename.c_str(), std::ofstream::app);
    outfile << t*simunits.t.outscale 																					<< "    ";
    outfile << Nsteps 																												<< "    ";
    for (k=0; k<ndim; k++) outfile << sink.star->r[k]*simunits.r.outscale 		<< "    ";
    for (k=0; k<ndim; k++) outfile << sink.star->v[k]*simunits.v.outscale 		<< "    ";
    for (k=0; k<ndim; k++) outfile << sink.star->a[k]*simunits.a.outscale 		<< "    ";
    for (k=0; k<ndim; k++) outfile << sink.angmom[k]*simunits.angmom.outscale << "    ";
    outfile << sink.star->m*simunits.m.outscale  															<< "    ";
    outfile << sink.menc     																									<< "    ";
    outfile << sink.mmax     																									<< "    ";
    outfile << sink.macctot  																									<< "    ";
    outfile << sink.dmdt*simunits.dmdt.outscale     													<< "    ";
    outfile << sink.ketot*simunits.E.outscale    															<< "    ";
    outfile << sink.gpetot*simunits.E.outscale   															<< "    ";
    outfile << sink.rotketot*simunits.E.outscale 															<< "    ";
    outfile << sink.utot*simunits.u.outscale     															<< "    ";
    outfile << sink.taccrete*simunits.t.outscale 															<< "    ";
    outfile << sink.trad*simunits.t.outscale     															<< "    ";
    outfile << sink.trot*simunits.t.outscale     															<< "    ";
    outfile << sink.tvisc*simunits.t.outscale    															<< "    ";
    for (k=0; k<ndim; k++) outfile << sink.L_star[k]*simunits.L.outscale    	<< "    ";
    outfile << sink.acclum*simunits.L.outscale    														<< "    ";
		outfile << sink.dmdt_star*simunits.dmdt.outscale   												<< "    ";
		outfile << sink.M_star*simunits.m.outscale    														<< "    ";
		outfile << sink.dmdt_MRI_0*simunits.dmdt.outscale													<< "    ";
		for (k=0; k<ndim; k++) outfile << sink.L_IAD[k]*simunits.angmom.outscale 	<< "    ";
		outfile << sink.stellarMd*simunits.m.outscale    													<< "    ";
		outfile << sink.stellarRadius*simunits.r.outscale    											<< "    ";
    outfile << sink.stellarPI    																							<< "    ";
    outfile << sink.stellarStage    																					<< "    ";
		outfile << sink.luminosity*simunits.L.outscale			   										<< "    ";
    outfile << sink.stellarTc*simunits.temp.outscale 													<< endl;
    outfile.close();

  }
  //-----------------------------------------------------------------------------------------------

  return;
}
