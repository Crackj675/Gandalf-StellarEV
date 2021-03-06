//=================================================================================================
//  Simulation.cpp
//  Contains all main functions controlling the simulation work-flow.
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
#include <sstream>
#include <string>
#include <math.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include "Precision.h"
#include "Exception.h"
#include "Debug.h"
#include "Simulation.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Nbody.h"
#include "Hydrodynamics.h"
#include "Sph.h"
#include "RiemannSolver.h"
#include "SimulationIO.hpp"
#include "SimulationIC.hpp"
#include "SimAnalysis.hpp"
#include "SphSnapshot.h"
using namespace std;


#ifdef _OPENMP
#include "omp.h"
#endif

// Declare invndim constant here (prevents warnings with some compilers)
template <int ndim>
const FLOAT Simulation<ndim>::invndim = 1.0/ndim;

#ifdef MPI_PARALLEL
const bool SimulationBase::MPI=true;
#else
const bool SimulationBase::MPI=false;
#endif

//=================================================================================================
//  SimulationBase::SimulationFactory
/// Creates a simulation object depending on the dimensionality.
//=================================================================================================
SimulationBase* SimulationBase::SimulationFactory
 (int ndim,                            ///< [in] No. of dimensions
  string simtype,                      ///< [in] Simulation type
  Parameters* params)                  ///< [in] Pointer to parameters object
{
  debug1("[SimulationBase::SimulationFactory]");

  // Check ndim is valid
  if (ndim < 1 || ndim > 3) {
    stringstream msg;
    msg << "Error: ndim must be either 1, 2, 3; the value " << ndim << "is not allowed!";
    ExceptionHandler::getIstance().raise(msg.str());
  }

  // Check simulation type is valid
  if (simtype != "sph" && simtype != "gradhsph" && simtype != "sm2012sph" &&
      simtype != "meshlessfv" && simtype != "mfvmuscl" && simtype != "mfvrk" &&
      simtype != "nbody" ) {
    string msg = "Error: the simulation type " + simtype + " was not recognized";
    ExceptionHandler::getIstance().raise(msg);
  }


  // Set ndim and simtype inside the parameters
  params->intparams["ndim"] = ndim;
  params->stringparams["sim"] = simtype;


  // Create and return Simulation object depending on the chosen algorithm
  // and the dimensionality.
  if (ndim == 1) {
    if (simtype == "gradhsph" || simtype == "sph") {
      return new GradhSphSimulation<1>(params);
    }
    else if (simtype == "sm2012sph") {
      return new SM2012SphSimulation<1>(params);
    }
    else if (simtype == "meshlessfv" || simtype == "mfvmuscl") {
      return new MfvMusclSimulation<1>(params);
    }
    else if (simtype == "mfvrk") {
      return new MfvRungeKuttaSimulation<1>(params);
    }
    else if (simtype == "nbody") {
      return new NbodySimulation<1>(params);
    }
  }
  else if (ndim == 2) {
    if (simtype == "gradhsph" || simtype == "sph") {
      return new GradhSphSimulation<2>(params);
    }
    else if (simtype == "sm2012sph") {
      return new SM2012SphSimulation<2>(params);
    }
    else if (simtype == "meshlessfv" || simtype == "mfvmuscl") {
      return new MfvMusclSimulation<2>(params);
    }
    else if (simtype == "mfvrk") {
      return new MfvRungeKuttaSimulation<2>(params);
    }
    else if (simtype == "nbody") {
      return new NbodySimulation<2>(params);
    }
  }
  else if (ndim == 3) {
    if (simtype == "gradhsph" || simtype == "sph") {
      return new GradhSphSimulation<3>(params);
    }
    else if (simtype == "sm2012sph") {
      return new SM2012SphSimulation<3>(params);
    }
    else if (simtype == "meshlessfv" || simtype == "mfvmuscl") {
      return new MfvMusclSimulation<3>(params);
    }
    else if (simtype == "mfvrk") {
      return new MfvRungeKuttaSimulation<3>(params);
    }
    else if (simtype == "nbody") {
      return new NbodySimulation<3>(params);
    }
  }
  return NULL;
}



//=================================================================================================
//  SimulationBase::SimulationBase
/// SimulationBase constructor, initialising important simulation variables.
//=================================================================================================
SimulationBase::SimulationBase
 (Parameters* params)                ///< [in] Pointer to parameters object
{
  simparams = new Parameters(*params);

  timing = NULL;
  level_max=0;
  ndims=1;
  dt_litesnap=0;
  tsnapfirst=0;
  recomputeRadiation=false;
  ntreestockstep=1;
  nbody_single_timestep=0;
  ntreebuildstep=1;
  dt_snap=0.2;
  noutputstep=128;
  ndiagstep=1024;
  tlitesnapnext=0;
  tend=1.0;
  level_step=1;
  hydro_single_timestep=0;
  level_diff_max=1;
  sink_particles=0;
  dt_max=0;
  nradstep=1;
  Nstepsmax=999999;
  extra_sink_output=false;
  dt_min_nbody=big_number_dp;
  dt_min_hydro=big_number_dp;
  pruning_level_max=6;
  rebuild_tree=false; //i think this was being used uninitialised
  tsnapnext=0;
  pruning_level_min=6;
  Nlevels=1;
  nsystembuildstep=1;
  dt_python=8.0;
  tmax_wallclock=1e20;


  paramfile             = "";
  integration_step      = 1;
  litesnap              = 0;
  n                     = 0;
  nlastrestart          = 0;
  nrestartstep          = 0;
  nresync               = 0;
  Nblocksteps           = 0;
  Nfullsteps            = 0;
  Nmpi                  = 1;
  Noutsnap              = 0;
  Noutlitesnap          = 0;
  Nsteps                = 0;
  rank                  = 0;
  dt_snap_wall          = 0.0;
  t                     = 0.0;
  timestep              = 0.0;
  tsnaplast             = 0.0;
  tlitesnaplast         = 0.0;
  tsnap_wallclock       = 0.0;
  ewaldGravity          = false;
  initial_h_provided    = false;
  kill_simulation       = false;
  ParametersProcessed   = false;
  periodicBoundaries    = false;
  rescale_particle_data = false;
  restart               = false;
  setup                 = false;
#if defined _OPENMP
  if (omp_get_dynamic()) {
    cout << "Warning: the dynamic adjustment of the number threads was on. "
         << "For better load-balancing, we will disable it" << endl;
  }
  omp_set_dynamic(0);
  Nthreads = omp_get_max_threads();
  assert(Nthreads > 0);
#else
  Nthreads = 1;
#endif
  timing = new CodeTiming() ;
}



//=================================================================================================
//  SimulationBase::~SimulationBase
/// SimulationBase destructor
//=================================================================================================
SimulationBase::~SimulationBase()
{
  if (timing != NULL)
    delete timing ;
}



//=================================================================================================
//  SimulationBase::SplashScreen
/// Write splash screen to standard output.
//=================================================================================================
void SimulationBase::SplashScreen(string& paramfile)
{
  cout << "******************************************************************************" << endl;
  cout << "*                                                                            *" << endl;
  cout << "*         *****     ****    *     *   *****     ****    *      ******        *" << endl;
  cout << "*        *     *   *    *   **    *   *    *   *    *   *      *             *" << endl;
  cout << "*        *         *    *   * *   *   *    *   *    *   *      *             *" << endl;
  cout << "*        *    **   ******   *  *  *   *    *   ******   *      ******        *" << endl;
  cout << "*        *     *   *    *   *   * *   *    *   *    *   *      *             *" << endl;
  cout << "*        *     *   *    *   *    **   *    *   *    *   *      *             *" << endl;
  cout << "*         *****    *    *   *     *   *****    *    *   *****  *             *" << endl;
  cout << "*                                                                            *" << endl;
  cout << "*   Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids    *" << endl;
  cout << "*                        Version 0.4.0 - 04/09/2015                          *" << endl;
  cout << "*                                                                            *" << endl;
  cout << "*                 Original code : D. A. Hubber & G. Rosotti                  *" << endl;
  cout << "*                                                                            *" << endl;
  cout << "*              Contributions by : S. Balfour, F. Dinnbier, S. Heigl,         *" << endl;
  cout << "*                                 O. Lomax, J. Ngoumou, P. Rohde,            *" << endl;
  cout << "*                                 S. Walch, A. P. Whitworth, R. Wunsch       *" << endl;
  cout << "*                                                                            *" << endl;
  cout << "*                  https://github.com/gandalfcode/gandalf                    *" << endl;
  cout << "*                                                                            *" << endl;
  cout << "******************************************************************************" << endl;
  cout << "Running from parameter file " << paramfile << endl;
  return;
}



//=================================================================================================
//  SimulationBase::SetParam
/// Accessor function for modifying a string value. Also checks that the
/// non return point has not been reached
//=================================================================================================
void SimulationBase::SetParam
 (string key,                          ///< [in] Parameter string name
  string value)                        ///< [in] Parameter string value
{
  // Error checking
  if (ParametersProcessed) {
    string msg = "Error: the non-return point for setting parameters has been reached!";
    ExceptionHandler::getIstance().raise(msg);
  }
  if (key == "ndim") {
    string msg = "Error: Not possible to change the number of dimensions!";
    ExceptionHandler::getIstance().raise(msg);
  }
  if (key == "sim") {
    string msg = "Error: Cannot change the type of simulation afterwards!";
    ExceptionHandler::getIstance().raise(msg);
  }

  simparams->SetParameter(key, value);
}



//=================================================================================================
//  SphSimulationBase::SetParam
/// Accessor function for modifying an int value, wrapper around the one for
/// string value. Also checks that the non return point has not been reached
//=================================================================================================
void SimulationBase::SetParam
 (string key,                          ///< [in] Parameter string name
  int value)                           ///< [in] Parameter integer value
{
  ostringstream convert;
  convert << value;
  SetParam (key, convert.str());
}



//=================================================================================================
//  SimulationBase::SetParam
/// Accessor function for modifying a float value, wrapper around the one for
/// string value.  Also checks that the non return point has not been reached.
//=================================================================================================
void SimulationBase::SetParam
 (string key,                          ///< [in] Parameter string name
  double value)                        ///< [in] Parameter float (double) value
{
  ostringstream convert;
  convert << value;
  SetParam (key, convert.str());
}



//=================================================================================================
//  SimulationBase::GetParam
/// Accessor function for getting a parameter value
/// Wrapper around the corresponding function in Parameters
//=================================================================================================
string SimulationBase::GetParam(string key)
{
  return simparams->GetParameter(key);
}



//=================================================================================================
//  SimulationBase::GetIntAndFloatParameterKeys
/// Returns a list containing the keys of all the int and float parameters
//=================================================================================================
std::list<string>* SimulationBase::GetIntAndFloatParameterKeys()
{
  if (! keys.empty()) return &keys;

  for (std::map<string, int>::iterator it=simparams->intparams.begin() ;
       it != simparams->intparams.end(); it++) {
    keys.push_back(it->first);
  }

  for (std::map<string, double>::iterator it=simparams->floatparams.begin() ;
       it != simparams->floatparams.end(); it++) {
    keys.push_back(it->first);
  }

  return &keys;
}



//=================================================================================================
//  SimulationBase::Run
/// Controls the simulation main loop, including exit conditions.  If provided as an optional
/// argument, will only advance the simulation by 'Nadvance' steps.
//=================================================================================================
void SimulationBase::Run
 (int Nadvance)                        ///< [in] Selected max no. of timesteps (Optional argument)
{
  int Ntarget;                         // Target step no before finishing main code integration.

  debug1("[SimulationBase::Run]");

  // Set integer timestep exit condition if provided as parameter.
  if (Nadvance < 0) Ntarget = Nstepsmax;
  else Ntarget = (int) (Nsteps + Nadvance);

  // Continue to run simulation until we reach the required time, or
  // exeeded the maximum allowed number of steps.
  //-----------------------------------------------------------------------------------------------
  while (t < tend && Nsteps < Ntarget) {

   CodeTiming::BlockTimer timer = timing->StartNewTimer("RUN");

    MainLoop();
    Output();

    // Special condition to check if maximum wall-clock time has been reached.
    if (kill_simulation || timing->RunningTime() > 0.95*tmax_wallclock) {
      RestartSnapshot();
      cout << "Reached maximum wall-clock time.  Killing simulation." << endl;
      break;
    }


  }
  //-----------------------------------------------------------------------------------------------

  FinaliseSimulation();
  CalculateDiagnostics();
  OutputDiagnostics();
  UpdateDiagnostics();

  cout << "Final t : " << t*simunits.t.outscale << " " << simunits.t.outunit
       << "    Total no. of steps : " << Nsteps << endl;


  // If reached end of simulation, remove 'cont' file to prevent automatic restart
  if (t >= tend || Nsteps >= Ntarget) {
    if (remove("cont") != 0) {
      cout << "Error deleting cont file" << endl;
    }
  }

  return;
}



//=================================================================================================
//  SimulationBase::InteractiveRun
/// Controls the simulation main loop, including exit conditions.
/// If provided, will only advance the simulation by 'Nadvance' steps.
//=================================================================================================
list<SphSnapshotBase*> SimulationBase::InteractiveRun
 (int Nadvance)                        ///< [in] Max no. of integer steps (Optional argument)
{
  int Ntarget;                // Selected integer timestep
  DOUBLE tdiff = 0.0;                  // Measured time difference
  clock_t tstart = clock();            // Initial CPU clock time
  string filename;                     // Name of the output file
  list<SphSnapshotBase*> snap_list;    // List of snapshots produced while running
                                       // that will be passed back to Python

  debug2("[SimulationBase::InteractiveRun]");

  // Set integer timestep exit condition if provided as parameter.
  if (Nadvance < 0) Ntarget = Nstepsmax;
  else Ntarget = (int) (Nsteps + Nadvance);

  // Continue to run simulation until we reach the required time, or
  // exeeded the maximum allowed number of steps.
  //-----------------------------------------------------------------------------------------------
  while (t < tend && Nsteps < Ntarget && tdiff < dt_python) {

    // Evolve the simulation one step
    MainLoop();

    // Update all diagnostics (including binaries) here for now
    if (t >= tsnapnext) CalculateDiagnostics();

    // Call output routine
    filename = Output();

    // If we have written a snapshot, create a new snapshot object
    if (filename.length() != 0) {
      SphSnapshotBase* snapshot =
        SphSnapshotBase::SphSnapshotFactory(filename, this, ndims);
      snapshot->CopyDataFromSimulation();
      snap_list.push_back(snapshot);
    }

    // Measure CPU clock time difference since current function was called
    tdiff = (DOUBLE) (clock() - tstart) / (DOUBLE) CLOCKS_PER_SEC;

  }
  //-----------------------------------------------------------------------------------------------


  // Calculate and process all diagnostic quantities
  if (t >= tend || Nsteps >= Ntarget) {
    FinaliseSimulation();
    CalculateDiagnostics();
    OutputDiagnostics();
    UpdateDiagnostics();
  }

  return snap_list;
}



//=================================================================================================
//  SimulationBase::Output
/// Controls when regular output snapshots are written by the code.
//=================================================================================================
string SimulationBase::Output(void)
{
  string filename;                  // 'Lite' output snapshot filename
  string filename2;                 // Regular output snapshot filename
  string nostring;                  // String of number of snapshots
  string fileend;                   // Name of restart file
  stringstream ss;                  // Stream object for preparing filename
  ofstream outfile;                 // Stream of restart file

  debug2("[SimulationBase::Output]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("OUTPUT");


  // Output time and no of steps for root process
  if (rank == 0) {
    if (Nsteps%noutputstep == 0) {
      cout << "t : " << t*simunits.t.outscale << " " << simunits.t.outunit
           << "    dt : " << timestep*simunits.t.outscale << " "
           << simunits.t.outunit << "    Nsteps : " << Nsteps << endl;
    }
  }


  // Output a lite-data snapshot for producing movies
  //-----------------------------------------------------------------------------------------------
  if (litesnap == 1 && t >= tlitesnapnext) {

    // Prepare filename for new snapshot
    Noutlitesnap++;
    tlitesnaplast = tlitesnapnext;
    tlitesnapnext += dt_litesnap;
    nostring = "";
    ss << setfill('0') << setw(5) << Noutlitesnap;
    nostring = ss.str();
    filename = run_id + ".slite." + nostring;
    ss.str(std::string());
    WriteSnapshotFile(filename,"slite");

  }
  //-----------------------------------------------------------------------------------------------

  // Output a data snapshot if reached required time
  //-----------------------------------------------------------------------------------------------
  if (t >= tsnapnext ) {

    // Prepare filename for new snapshot
    Noutsnap++;
    tsnaplast = tsnapnext;
    tsnapnext += dt_snap;
    nostring = "";
    ss << setfill('0') << setw(5) << Noutsnap;
    nostring = ss.str();
    filename = run_id + '.' + out_file_form + '.' + nostring;
    ss.str(std::string());
    WriteSnapshotFile(filename,out_file_form);

    // Now write name and format of snapshot to file (for restarts)
    if (rank == 0) {

      fileend = "restart";
      filename2 = run_id + "." + fileend;
      outfile.open(filename2.c_str());
      outfile << out_file_form << endl;
      outfile << filename << endl;
      outfile.close();

      // Finally, calculate wall-clock time interval since last output snapshot
      if (tsnap_wallclock > 0.0) dt_snap_wall = timing->RunningTime() - tsnap_wallclock;
      tsnap_wallclock = timing->RunningTime();

      // If simulation is too close to maximum wall-clock time, end prematurely
      if (timing->RunningTime() > 0.95*tmax_wallclock) {
        kill_simulation = true;
      }

    }

  }
  //-----------------------------------------------------------------------------------------------


  // Output diagnostics to screen if passed sufficient number of block steps
  if (Nblocksteps%ndiagstep == 0 && n%nresync == 0) {
    CalculateDiagnostics();
    OutputDiagnostics();
    UpdateDiagnostics();
    timing->ComputeTimingStatistics(run_id);

  }

  // Create temporary snapshot file
  if (n%nresync == 0 && Nsteps - nlastrestart >= nrestartstep) {
    RestartSnapshot();
    nlastrestart = Nsteps;
  }


  return filename;
}



//=================================================================================================
//  SimulationBase::RestartSnapshot
/// Write the restart log file (containing the last snapshot i.d.) plus a temporary snapshot
/// file for future restarting.
//=================================================================================================
void SimulationBase::RestartSnapshot(void)
{
  string filename;                     // Temporary output snapshot filename
  string filename2;                    // Restart log filename
  stringstream ss;                     // Stream object for preparing filename
  ofstream outfile;                    // Stream of restart file

  debug2("[SimulationBase::RestartSnapshot]");

  // Prepare filename for new snapshot
  filename = run_id + "." + out_file_form + ".tmp";
  ss.str(std::string());
  WriteSnapshotFile(filename,out_file_form);

  // Now write name and format of snapshot to file (for restarts)
  filename2 = run_id + ".restart";
  outfile.open(filename2.c_str());
  outfile << out_file_form << endl;
  outfile << filename << endl;
  outfile.close();

  return;
}



//=================================================================================================
//  SimulationBase::SetupSimulation
/// Main function for setting up a new simulation.
//=================================================================================================
void SimulationBase::SetupSimulation(void)
{
  debug1("[SimulationBase::Setup]");

  CodeTiming::BlockTimer timer = timing->StartNewTimer("SETUP");

  if (setup) {
    string msg = "This simulation has been already set up";
    ExceptionHandler::getIstance().raise(msg);
  }

  // Process the parameters file setting up all simulation objects
  if (simparams->stringparams["ic"] == "python") {
    if (!ParametersProcessed) {
      string msg = "Error: you are attempting to setup a simulation with initial conditions "
                   "generated from Python. Before setting up the simulation, you need to "
                   "import the initial conditions";
      ExceptionHandler::getIstance().raise(msg);
    }
  }
  else {
    if (ParametersProcessed) {
      string msg = "The parameters of the simulation have been already processed. It means that "
                   "you shouldn't be calling this function, please consult the documentation.";
      ExceptionHandler::getIstance().raise(msg);
    }
    ProcessParameters();
  }

  // Generate initial conditions for simulation on root process (for MPI jobs)
  if (rank == 0) {
    GenerateIC();
  }

  // Change to COM frame if selected
  if (simparams->intparams["com_frame"] == 1) SetComFrame();

  // Perform the rest of the initialisation, calculating all initial particle
  // quantities and setting up trees.
  PostInitialConditionsSetup();

  // Initial output before simulation begins
  Output();

  return;
}



//=================================================================================================
//  Simulation::ProcessNbodyParameters
/// Process all the options chosen in the parameters file for setting up
/// objects related to stars and N-body integration.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::ProcessNbodyParameters(void)
{
  map<string, int> &intparams = simparams->intparams;
  map<string, double> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;

  // Create N-body object based on chosen method in params file
  //-----------------------------------------------------------------------------------------------
  if (stringparams["nbody"] == "lfkdk") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyLeapfrogKDK<ndim, TabulatedKernel>
        (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
         floatparams["nbody_mult"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      if (KernelName == "m4") {
        nbody = new NbodyLeapfrogKDK<ndim, M4Kernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "quintic") {
        nbody = new NbodyLeapfrogKDK<ndim, QuinticKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "gaussian") {
        nbody = new NbodyLeapfrogKDK<ndim, GaussianKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else {
        string message = "Unrecognised parameter : kernel = " +
          simparams->stringparams["kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
        stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  //-----------------------------------------------------------------------------------------------
  else if (stringparams["nbody"] == "lfdkd") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyLeapfrogDKD<ndim, TabulatedKernel>
        (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
         floatparams["nbody_mult"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      if (KernelName == "m4") {
        nbody = new NbodyLeapfrogDKD<ndim, M4Kernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "quintic") {
        nbody = new NbodyLeapfrogDKD<ndim, QuinticKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "gaussian") {
        nbody = new NbodyLeapfrogDKD<ndim, GaussianKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else {
        string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
        stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
    integration_step = max(integration_step,2);
  }
  //-----------------------------------------------------------------------------------------------
  else if (stringparams["nbody"] == "hermite4") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyHermite4<ndim, TabulatedKernel>
        (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
         floatparams["nbody_mult"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      if (KernelName == "m4") {
        nbody = new NbodyHermite4<ndim, M4Kernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "quintic") {
        nbody = new NbodyHermite4<ndim, QuinticKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "gaussian") {
        nbody = new NbodyHermite4<ndim, GaussianKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else {
        string message = "Unrecognised parameter : kernel = " +
          simparams->stringparams["kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
        stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  //-----------------------------------------------------------------------------------------------
  else if (stringparams["nbody"] == "hermite4ts") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyHermite4TS<ndim, TabulatedKernel>
        (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
         floatparams["nbody_mult"], KernelName, intparams["Npec"]);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      if (KernelName == "m4") {
        nbody = new NbodyHermite4TS<ndim, M4Kernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else if (KernelName == "quintic") {
        nbody = new NbodyHermite4TS<ndim, QuinticKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else if (KernelName == "gaussian") {
        nbody = new NbodyHermite4TS<ndim, GaussianKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else {
        string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
        stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  //-----------------------------------------------------------------------------------------------
  else if (stringparams["nbody"] == "hermite6ts") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyHermite6TS<ndim, TabulatedKernel>
        (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
         floatparams["nbody_mult"], KernelName, intparams["Npec"]);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      if (KernelName == "m4") {
        nbody = new NbodyHermite6TS<ndim, M4Kernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else if (KernelName == "quintic") {
        nbody = new NbodyHermite6TS<ndim, QuinticKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else if (KernelName == "gaussian") {
        nbody = new NbodyHermite6TS<ndim, GaussianKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else {
        string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
        stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  //-----------------------------------------------------------------------------------------------
  else {
    string message = "Unrecognised parameter : nbody = "
      + simparams->stringparams["nbody"];
    ExceptionHandler::getIstance().raise(message);
  }
  //-----------------------------------------------------------------------------------------------


  // Create sub-system object based on chosen method in params file
  //-----------------------------------------------------------------------------------------------
  if (intparams["sub_systems"] == 1) {

    //---------------------------------------------------------------------------------------------
    if (stringparams["sub_system_integration"] == "lfkdk") {
      string KernelName = stringparams["kernel"];
      if (intparams["tabulated_kernel"] == 1) {
        subsystem = new NbodyLeapfrogKDK<ndim, TabulatedKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["subsys_mult"], KernelName);
      }
      else if (intparams["tabulated_kernel"] == 0) {
        if (KernelName == "m4") {
          subsystem = new NbodyLeapfrogKDK<ndim, M4Kernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName);
        }
        else if (KernelName == "quintic") {
          subsystem = new NbodyLeapfrogKDK<ndim, QuinticKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName);
        }
        else if (KernelName == "gaussian") {
          subsystem = new NbodyLeapfrogKDK<ndim, GaussianKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName);
        }
        else {
          string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
          ExceptionHandler::getIstance().raise(message);
        }
      }
      else {
        string message = "Invalid option for the tabulated_kernel parameter: "
          + stringparams["tabulated_kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    //---------------------------------------------------------------------------------------------
    else if (stringparams["sub_system_integration"] == "hermite4") {
      string KernelName = stringparams["kernel"];
      if (intparams["tabulated_kernel"] == 1) {
        subsystem = new NbodyHermite4<ndim, TabulatedKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["subsys_mult"], KernelName);
      }
      else if (intparams["tabulated_kernel"] == 0) {
        if (KernelName == "m4") {
          subsystem = new NbodyHermite4<ndim, M4Kernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName);
        }
        else if (KernelName == "quintic") {
          subsystem = new NbodyHermite4<ndim, QuinticKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName);
        }
        else if (KernelName == "gaussian") {
          subsystem = new NbodyHermite4<ndim, GaussianKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName);
        }
        else {
          string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
          ExceptionHandler::getIstance().raise(message);
        }
      }
      else {
        string message = "Invalid option for the tabulated_kernel parameter: "
          + stringparams["tabulated_kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    //---------------------------------------------------------------------------------------------
    else if (stringparams["sub_system_integration"] == "hermite4ts") {
      string KernelName = stringparams["kernel"];
      if (intparams["tabulated_kernel"] == 1) {
        subsystem = new NbodyHermite4TS<ndim, TabulatedKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["subsys_mult"], KernelName, intparams["Npec"]);
      }
      else if (intparams["tabulated_kernel"] == 0) {
        if (KernelName == "m4") {
          subsystem = new NbodyHermite4TS<ndim, M4Kernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName, intparams["Npec"]);
        }
        else if (KernelName == "quintic") {
          subsystem = new NbodyHermite4TS<ndim, QuinticKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName, intparams["Npec"]);
        }
        else if (KernelName == "gaussian") {
          subsystem = new NbodyHermite4TS<ndim, GaussianKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName, intparams["Npec"]);
        }
        else {
          string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
          ExceptionHandler::getIstance().raise(message);
        }
      }
      else {
        string message = "Invalid option for the tabulated_kernel parameter: "
          + stringparams["tabulated_kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    //---------------------------------------------------------------------------------------------
    else if (stringparams["sub_system_integration"] == "hermite6ts") {
      string KernelName = stringparams["kernel"];
      if (intparams["tabulated_kernel"] == 1) {
        subsystem = new NbodyHermite6TS<ndim, TabulatedKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["subsys_mult"], KernelName, intparams["Npec"]);
      }
      else if (intparams["tabulated_kernel"] == 0) {
        if (KernelName == "m4") {
          subsystem = new NbodyHermite6TS<ndim, M4Kernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName, intparams["Npec"]);
        }
        else if (KernelName == "quintic") {
          subsystem = new NbodyHermite6TS<ndim, QuinticKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName, intparams["Npec"]);
        }
        else if (KernelName == "gaussian") {
          subsystem = new NbodyHermite6TS<ndim, GaussianKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName, intparams["Npec"]);
        }
        else {
          string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
          ExceptionHandler::getIstance().raise(message);
        }
      }
      else {
        string message = "Invalid option for the tabulated_kernel parameter: "
          + stringparams["tabulated_kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    //---------------------------------------------------------------------------------------------
    else {
      string message = "Unrecognised parameter : sub_system_integration = "
        + simparams->stringparams["sub_system_integration"];
      ExceptionHandler::getIstance().raise(message);
    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}

//=================================================================================================
//  Simulation::ProcessParameters
/// Process all the options chosen in the parameters file, setting various
/// simulation variables and creating important simulation objects.
/// This initialises all the things in common between different schemes
//=================================================================================================
template <int ndim>
void Simulation<ndim>::ProcessParameters(void)
{

  debug2("[Simulation::ProcessParameters]");

  map<string, string> &stringparams = simparams->stringparams;
  map<string, int> &intparams = simparams->intparams;
  map<string, double> &floatparams = simparams->floatparams;


  // Sanity check for valid dimensionality
  if (ndim < 1 || ndim > 3) {
    std::ostringstream message;
    message << "Invalid dimensionality chosen : ndim = " << ndim;
    ExceptionHandler::getIstance().raise(message.str());
  }
  // Set-up random number generator object
  //-----------------------------------------------------------------------------------------------
  if (stringparams["rand_algorithm"] == "xorshift") {
    randnumb = new XorshiftRand(intparams["randseed"]);
  }
  else if (stringparams["rand_algorithm"] == "none") {
    randnumb = new DefaultSystemRand(intparams["randseed"]);
  }
  else {
    string message = "Unrecognised parameter : rand_algorithm= " +
      stringparams["rand_algorithm"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Set-up all output units for scaling parameters
  simunits.SetupUnits(simparams);

  // Boundary condition variables
  //-----------------------------------------------------------------------------------------------
  simbox.boundary_lhs[0] = setBoundaryType(stringparams["boundary_lhs[0]"]);
  simbox.boundary_rhs[0] = setBoundaryType(stringparams["boundary_rhs[0]"]);
  if (simbox.boundary_lhs[0] != openBoundary) simbox.min[0] = floatparams["boxmin[0]"]/simunits.r.outscale;
  if (simbox.boundary_rhs[0] != openBoundary) simbox.max[0] = floatparams["boxmax[0]"]/simunits.r.outscale;

  simbox.boundary_lhs[1] = setBoundaryType(stringparams["boundary_lhs[1]"]);
  simbox.boundary_rhs[1] = setBoundaryType(stringparams["boundary_rhs[1]"]);
  if (simbox.boundary_lhs[1] != openBoundary) simbox.min[1] = floatparams["boxmin[1]"]/simunits.r.outscale;
  if (simbox.boundary_rhs[1] != openBoundary) simbox.max[1] = floatparams["boxmax[1]"]/simunits.r.outscale;

  simbox.boundary_lhs[2] = setBoundaryType(stringparams["boundary_lhs[2]"]);
  simbox.boundary_rhs[2] = setBoundaryType(stringparams["boundary_rhs[2]"]);
  if (simbox.boundary_lhs[2] != openBoundary) simbox.min[2] = floatparams["boxmin[2]"]/simunits.r.outscale;
  if (simbox.boundary_rhs[2] != openBoundary) simbox.max[2] = floatparams["boxmax[2]"]/simunits.r.outscale;

  for (int k=0; k<ndim; k++) {
    simbox.size[k] = simbox.max[k] - simbox.min[k];
    simbox.half[k] = (FLOAT) 0.5*simbox.size[k];
  }


  // Initial conditions box variables
  //-----------------------------------------------------------------------------------------------
  for (int k=0; k<ndim; k++) icBox.boundary_lhs[k] = periodicBoundary;
  for (int k=0; k<ndim; k++) icBox.boundary_rhs[k] = periodicBoundary;
  icBox.min[0] = simparams->floatparams["boxmin[0]"]/simunits.r.outscale;
  icBox.max[0] = simparams->floatparams["boxmax[0]"]/simunits.r.outscale;
  icBox.min[1] = simparams->floatparams["boxmin[1]"]/simunits.r.outscale;
  icBox.max[1] = simparams->floatparams["boxmax[1]"]/simunits.r.outscale;
  icBox.min[2] = simparams->floatparams["boxmin[2]"]/simunits.r.outscale;
  icBox.max[2] = simparams->floatparams["boxmax[2]"]/simunits.r.outscale;
  for (int k=0; k<ndim; k++) {
    icBox.size[k] = icBox.max[k] - icBox.min[k];
    icBox.half[k] = (FLOAT) 0.5*icBox.size[k];
  }



  // Set external potential field object
  if (stringparams["external_potential"] == "none") {
    extpot = new NullPotential<ndim>();
  }
  else if (stringparams["external_potential"] == "vertical") {
    extpot = new VerticalPotential<ndim>
      (intparams["kgrav"], floatparams["avert"], simbox.min[intparams["kgrav"]]);
  }
  else if (stringparams["external_potential"] == "plummer") {
    extpot = new PlummerPotential<ndim>(floatparams["mplummer"], floatparams["rplummer"]);
  }
  else if (stringparams["external_potential"] == "silcc") {
    extpot = new SilccPotential<ndim>(floatparams["sigma_star"], floatparams["z_d"], simunits);
  }
  else {
    string message = "Unrecognised parameter : external_potential = "
      + simparams->stringparams["external_potential"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Create Ewald periodic gravity object
  periodicBoundaries = IsAnyBoundaryPeriodic(simbox);
  if (IsAnyBoundaryReflecting(simbox) && intparams["self_gravity"]) {
    ExceptionHandler::getIstance().raise("Error: Reflecting boundaries and self-gravity is not "
                                         "supported") ;
  }
  else if (periodicBoundaries && intparams["self_gravity"] == 1) {
    if (ndim != 3) {
      ExceptionHandler::getIstance().raise("Error: Periodic/Ewald gravity only supported in 3D");
    }
    ewaldGravity = true;
    ewald = new Ewald<ndim>
      (simbox, intparams["gr_bhewaldseriesn"], intparams["in"], intparams["nEwaldGrid"],
       floatparams["ewald_mult"], floatparams["ixmin"], floatparams["ixmax"],
       floatparams["EFratio"], timing);
    simbox.PeriodicGravity = true ;
  }
  else {
    simbox.PeriodicGravity = false ;
  }

  // Set other important simulation variables
  dt_litesnap         = floatparams["dt_litesnap"]/simunits.t.outscale;
  dt_python           = floatparams["dt_python"];
  dt_snap             = floatparams["dt_snap"]/simunits.t.outscale;
  extra_sink_output   = intparams["extra_sink_output"];
  level_diff_max      = intparams["level_diff_max"];
  litesnap            = intparams["litesnap"];
  Nlevels             = intparams["Nlevels"];
  ndiagstep           = intparams["ndiagstep"];
  noutputstep         = intparams["noutputstep"];
  nrestartstep        = intparams["nrestartstep"];
  ntreebuildstep      = intparams["ntreebuildstep"];
  ntreestockstep      = intparams["ntreestockstep"];
  nsystembuildstep    = intparams["nsystembuildstep"];
  Nstepsmax           = intparams["Nstepsmax"];
  out_file_form       = stringparams["out_file_form"];
  pruning_level_min   = intparams["pruning_level_min"];
  pruning_level_max   = intparams["pruning_level_max"];
  run_id              = stringparams["run_id"];
  hydro_single_timestep = intparams["sph_single_timestep"];
  tmax_wallclock      = floatparams["tmax_wallclock"];
  tend                = floatparams["tend"]/simunits.t.outscale;
  tlitesnapnext       = floatparams["tlitesnapfirst"]/simunits.t.outscale;
  tsnapnext           = floatparams["tsnapfirst"]/simunits.t.outscale;

}

//=================================================================================================
//  Simulation::AllocateParticleMemory
/// Allocate all memory for both SPH and N-body particles.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::AllocateParticleMemory(void)
{
  int N;                               // Max. no. of stars/sinks

  debug2("[Simulation::AllocateParticleMemory]");

  // Allocate N-body memory (if using N-body)
  //-----------------------------------------------------------------------------------------------
  if (nbody) {

    // If sink particles are employed, allow enough memory for new sinks
    if (sink_particles == 1) {
      N = max(nbody->Nstar, 1024);
    }
    else N = nbody->Nstar;

    // Now call all memory allocation routines
    nbody->AllocateMemory(N);
    sinks->AllocateMemory(N);
  }
  //-----------------------------------------------------------------------------------------------

  // Allocate SPH memory, if being used
  if (hydro) hydro->AllocateMemory(hydro->Nhydro);


  return;
}



//=================================================================================================
//  Simulation::DeallocateParticleMemory
/// Deallocate all particle memory
//=================================================================================================
template <int ndim>
void Simulation<ndim>::DeallocateParticleMemory(void)
{
  debug2("[Simulation::DellocateParticleMemory]");

  sinks->DeallocateMemory();
  nbody->DeallocateMemory();
  hydro->DeallocateMemory();

  return;
}



//=================================================================================================
//  Simulation::PreSetupForPython
/// Initialisation routine called by python interface.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::PreSetupForPython(void)
{
  debug1("[Simulation::PreSetupForPython]");

  // Check that IC type is really python
  if (simparams->stringparams["ic"] != "python") {
    string msg = "Error: you should call this function only if you are "
      "using \"python\" as \"ic\" parameter";
    ExceptionHandler::getIstance().raise(msg);
  }

  if (ParametersProcessed) {
    string msg = "Error: ProcessParameters has been already called!";
    ExceptionHandler::getIstance().raise(msg);
  }

  // Parse all parameters and set-up all objects required for simulation
  ProcessParameters();

  // Allocate all memory for both hydro and N-body particles
  hydro->Nhydro = simparams->intparams["Nhydro"];
  if (nbody)
    nbody->Nstar = simparams->intparams["Nstar"];
  AllocateParticleMemory();

  return;
}



//=================================================================================================
//  Simulation::ImportArrayNbody
/// Import an array containing nbody particle properties from python to C++ arrays.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::ImportArrayNbody
 (double* input,                       ///< [in] Array of values from python
  int size,                            ///< [in] No. of array elements
  string quantity)                     ///< [in] String id of quantity being imported
{
  FLOAT StarParticle<ndim>::*quantityp = 0;            // Pointer to scalar quantity
  FLOAT (StarParticle<ndim>::*quantitypvec)[ndim] = 0; // Pointer to component of vector quantity
  int index = 0;                                       // Component index (if quantity is vector)
  bool scalar = false;                                 // Is the requested quantity a scalar?

  // Check that the size is correct
  if (size != nbody->Nstar) {
    stringstream message;
    message << "Error: the array you are passing has a size of " << size
            << ", but memory has been allocated for " << nbody->Nstar << " star particles";
    ExceptionHandler::getIstance().raise(message.str());
  }

  // Now set pointer to the correct value inside the particle data structure
  //-----------------------------------------------------------------------------------------------
  if (quantity == "x") {
    quantitypvec = &StarParticle<ndim>::r;
    index = 0;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "y") {
    if (ndim < 2) {
      string message = "Error: loading y-coordinate array for ndim < 2";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &StarParticle<ndim>::r;
    index = 1;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "z") {
    if (ndim < 3) {
      string message = "Error: loading y-coordinate array for ndim < 3";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &StarParticle<ndim>::r;
    index = 2;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "vx") {
    quantitypvec = &StarParticle<ndim>::v;
    index = 0;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "vy") {
    if (ndim < 2) {
      string message = "Error: loading vy array for ndim < 2";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &StarParticle<ndim>::v;
    index = 1;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "vz") {
    if (ndim < 3) {
      string message = "Error: loading vz array for ndim < 3";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &StarParticle<ndim>::v;
    index = 2;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "m") {
    quantityp = &StarParticle<ndim>::m;
    scalar = true;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "h") {
    quantityp = &StarParticle<ndim>::h;
    scalar = true;
  }
  //-----------------------------------------------------------------------------------------------
  else {
    string message = "Quantity " + quantity + "not recognised";
    ExceptionHandler::getIstance().raise(message);
  }
  //-----------------------------------------------------------------------------------------------


  // Finally loop over particles and set all values
  // (Note that the syntax for scalar is different from the one for vectors)
  //-----------------------------------------------------------------------------------------------
  if (scalar) {
    int i=0;
    for (StarParticle<ndim>* particlep = nbody->stardata;
         particlep < nbody->stardata+size; particlep++, i++) {
      particlep->*quantityp = input[i];
    }
  }
  else {
    int i=0;
    for (StarParticle<ndim>* particlep = nbody->stardata;
         particlep < nbody->stardata+size; particlep++, i++) {
      (particlep->*quantitypvec)[index] = input[i];
    }
  }

  return;
}



//=================================================================================================
//  Simulation::ImportArraySph
/// Import an array containing sph particle properties from python to C++ arrays.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::ImportArraySph
 (double* input,                       ///< [in] Array of values imported from python
  int size,                            ///< [in] No. of elements in array
  string quantity)                     ///< [in] String id of quantity
{
  FLOAT Particle<ndim>::*quantityp = 0;             // Pointer to scalar quantity
  FLOAT (Particle<ndim>::*quantitypvec)[ndim] = 0;  // Pointer to component of vector quantity
  int index = 0;                                    // If it's a component of a vector
                                                    // quantity, we need to know its index
  bool scalar = false;                              // Is the requested quantity a scalar?

  // Check that the size is correct
  if (size != hydro->Nhydro) {
    stringstream message;
    message << "Error: the array you are passing has a size of "
            << size << ", but memory has been allocated for " << hydro->Nhydro << " particles";
    ExceptionHandler::getIstance().raise(message.str());
  }


  // Now set pointer to the correct value inside the particle data structure
  //-----------------------------------------------------------------------------------------------
  if (quantity == "x") {
    quantitypvec = &Particle<ndim>::r;
    index = 0;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "y") {
    if (ndim < 2) {
      string message = "Error: loading y-coordinate array for ndim < 2";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &Particle<ndim>::r;
    index = 1;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "z") {
    if (ndim < 3) {
      string message = "Error: loading y-coordinate array for ndim < 3";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &Particle<ndim>::r;
    index = 2;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "vx") {
    quantitypvec = &Particle<ndim>::v;
    index = 0;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "vy") {
    if (ndim < 2) {
      string message = "Error: loading vy array for ndim < 2";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &Particle<ndim>::v;
    index = 1;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "vz") {
    if (ndim < 3) {
      string message = "Error: loading vz array for ndim < 3";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &Particle<ndim>::v;
    index = 2;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "rho") {
    //TODO: at the moment, if rho or h are uploaded, they will be just ignored.
    //Add some facility to use them
    quantityp = &Particle<ndim>::rho;
    scalar = true;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "h") {
    quantityp = &Particle<ndim>::h;
    scalar = true;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "u") {
    //TODO: add some facility for uploading either u, T, or cs, and compute automatically the
    //other ones depending on the EOS
    quantityp = &Particle<ndim>::u;
    scalar=true;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "m") {
    quantityp = &Particle<ndim>::m;
    scalar = true;
  }
  //-----------------------------------------------------------------------------------------------
  else {
    string message = "Quantity " + quantity + "not recognised";
    ExceptionHandler::getIstance().raise(message);
  }
  //-----------------------------------------------------------------------------------------------


  // Finally loop over particles and set all values
  // (Note that the syntax for scalar is different from the one for vectors)
  //-----------------------------------------------------------------------------------------------
  if (scalar) {
    for (int i=0; i<size; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      part.*quantityp = input[i];
    }
  }
  else {
    for (int i=0; i<size; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      (part.*quantitypvec)[index] = input[i];
    }
  }

  return;
}



//=================================================================================================
//  Simulation::ImportArray
/// Import an array containing particle properties from python to C++ arrays.
/// This is a wrapper around ImportArraySph and ImportArrayNbody
//=================================================================================================
template <int ndim>
void Simulation<ndim>::ImportArray
 (double* input,                       ///< [in] Input array
  int size,                            ///< [in] Size of the input array
  string quantity,                     ///< [in] Quantity to be set equal to the given array
  string type)                         ///< [in] Particle type that should be assigned the array
{
  debug2("[Simulation::ImportArray]");

  // Check that PreSetup has been called
  if (! ParametersProcessed) {
    string msg = "Error: before calling ImportArray, you need to call PreSetupForPython!";
    ExceptionHandler::getIstance().raise(msg);
  }

  // Call the right function depending on the passed in type
  if (type == "sph") {
    // Check sph has been allocated
    if (hydro == NULL) {
      string message = "Error: memory for sph was not allocated! Are you sure that this is not a nbody-only simulation?";
      ExceptionHandler::getIstance().raise(message);
    }
    ImportArraySph(input, size, quantity);

  }
  else if (type == "star") {
    if (nbody == NULL) {
      string message = "Error: memory for nbody was not allocated! Are you sure that this is not a sph-only simulation?";
      ExceptionHandler::getIstance().raise(message);
    }
    ImportArrayNbody(input, size, quantity);
  }
  else {
    string message = "Error: we did not recognize the type " + type +
      ", the only allowed types are \"sph\" and \"nbody\"";
    ExceptionHandler::getIstance().raise(message);
  }

  return;
}



//=================================================================================================
//  Simulation::SetComFrame
/// Move all particles (both hydro and N-body) to centre-of-mass frame.
//=================================================================================================
template<int ndim>
void Simulation<ndim>::SetComFrame(void)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter

  debug2("[Simulation::SetComFrame]");

  CalculateDiagnostics();

  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.r[k] -= diag.rcom[k];
    for (k=0; k<ndim; k++) part.v[k] -= diag.vcom[k];
  }

  for (i=0; i<nbody->Nstar; i++) {
    for (k=0; k<ndim; k++) nbody->stardata[i].r[k] -= diag.rcom[k];
    for (k=0; k<ndim; k++) nbody->stardata[i].v[k] -= diag.vcom[k];
  }

  CalculateDiagnostics();

  return;
}



//=================================================================================================
//  Simulation::UpdateDiagnostics
/// Update energy error value after computing diagnostic quantities.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::UpdateDiagnostics(void)
{
  if (rank == 0) {
    diag.Eerror = fabs(diag0.Etot - diag.Etot)/fabs(diag0.Etot);
    cout << "Eerror : " << diag.Eerror << endl;
  }
}



//=================================================================================================
//  Simulation::ComputeGlobalTimestep
/// Computes global timestep for the simulation.  Calculates the minimum
/// timestep for all hydro and N-body particles in the simulation.
//=================================================================================================
template<int ndim>
void Simulation<ndim>::ComputeGlobalTimestep() {
  int i;                               // Particle counter
  DOUBLE dt;                           // Particle timestep
  DOUBLE dt_min = big_number_dp;       // Local copy of minimum timestep
  DOUBLE dt_nbody;                     // Aux. minimum N-body timestep
  DOUBLE dt_hydro;                       // Aux. minimum hydro timestep

  debug2("[Simulation::DoComputeGlobalTimestep]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("GLOBAL_TIMESTEPS");


   // Only update timestep when all particles are synced at end of last step.
   //-----------------------------------------------------------------------------------------------
   if (n == nresync) {

     n            = 0;
     level_max    = 0;
     level_step   = level_max + integration_step - 1;
     nresync      = integration_step;
     dt_min_nbody = big_number_dp;
     dt_min_hydro = big_number_dp;

     // Find minimum timestep from all hydro particles
     //---------------------------------------------------------------------------------------------
 #pragma omp parallel default(none) private(i,dt,dt_nbody,dt_hydro) shared(dt_min)
     {
       dt       = big_number_dp;
       dt_nbody = big_number_dp;
       dt_hydro = big_number_dp;

       if (hydro != NULL) {
 #pragma omp for
         for (i=0; i<hydro->Nhydro; i++) {
           Particle<ndim>& part = hydro->GetParticlePointer(i);
           part.flags.set(end_timestep) ;
           part.level     = 0;
           part.levelneib = 0;
           part.nstep     = pow(2,level_step - part.level);
           part.dt_next   = hydroint->Timestep(part,hydro);
           dt             = min(dt, (DOUBLE) part.dt_next);
           dt_hydro       = min(dt_hydro, (DOUBLE) part.dt_next);
         }
       }

       // Now compute minimum timestep due to stars/systems
       if (nbody != NULL) {
 #pragma omp for
         for (i=0; i<nbody->Nnbody; i++) {
           nbody->nbodydata[i]->flags.set(end_timestep);
           nbody->nbodydata[i]->level = 0;
           nbody->nbodydata[i]->nstep = pow(2,level_step - nbody->nbodydata[i]->level);
           nbody->nbodydata[i]->dt_next  = nbody->Timestep(nbody->nbodydata[i]);
           dt       = min(dt,nbody->nbodydata[i]->dt_next);
           dt_nbody = min(dt_nbody,nbody->nbodydata[i]->dt_next);
         }
       }

 #pragma omp critical
       {
         if (dt < dt_min) dt_min = dt;
         if (dt_hydro < dt_min_hydro) dt_min_hydro = dt_hydro;
         if (dt_nbody < dt_min_nbody) dt_min_nbody = dt_nbody;
       }

     }
     //---------------------------------------------------------------------------------------------

 #ifdef MPI_PARALLEL
     dt = dt_min;
     MPI_Allreduce(&dt, &dt_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
 #endif
     timestep = dt_min;

     // Set minimum timestep for all hydro and N-body particles
     if (hydro != NULL)
       for (i=0; i<hydro->Nhydro; i++) hydro->GetParticlePointer(i).dt_next = timestep;
     if (nbody != NULL)
       for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->dt_next = timestep;

   }
   //-----------------------------------------------------------------------------------------------

   return;


}




//=================================================================================================
//  Simulation::ComputeBlockTimesteps
/// Compute timesteps for all particles using hierarchical block timesteps.
//=================================================================================================
template<int ndim>
void Simulation<ndim>::ComputeBlockTimesteps(void)
{
  int i;                                     // Particle counter
  int istep;                                 // Aux. variable for changing steps
  int last_level;                            // Previous timestep level
  int level;                                 // Particle timestep level
  int level_max_aux;                         // Aux. maximum level variable
  int level_max_nbody = 0;                   // level_max for star particles only
  int level_max_old;                         // Old level_max
  int level_max_hydro = 0;                     // level_max for hydro particles only
  int level_nbody;                           // local thread var. for N-body level
  int level_hydro;                             // local thread var. for hydro level
  int nfactor;                               // Increase/decrease factor of n
  int nstep;                                 // Particle integer step-size
  DOUBLE dt;                                 // Aux. timestep variable
  DOUBLE dt_min = big_number_dp;             // Minimum timestep
  DOUBLE dt_min_aux;                         // Aux. minimum timestep variable
  DOUBLE dt_nbody;                           // Aux. minimum N-body timestep
  DOUBLE dt_hydro;                             // Aux. minimum hydro timestep
  DOUBLE timestep_temp[Nthreads];
  DOUBLE dt_min_hydro_temp[Nthreads];
  DOUBLE dt_min_nbody_temp[Nthreads];
  int level_max_temp[Nthreads];


  debug2("[Simulation::DoComputeBlockTimesteps]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("BLOCK_TIMESTEPS");


  dt_min_nbody = big_number_dp;
  dt_min_hydro = big_number_dp;



  // Synchronise all timesteps and reconstruct block timestep structure.
  //===============================================================================================
  if (n == nresync) {

    n = 0;
    timestep = big_number_dp;


#pragma omp parallel default(none) private(dt,dt_min_aux,dt_nbody,dt_hydro,i)\
  shared(timestep_temp,dt_min_hydro_temp,dt_min_nbody_temp)
    {
      // Initialise all timestep and min/max variables
      dt_min_aux = big_number_dp;
      dt_hydro   = big_number_dp;
      dt_nbody   = big_number_dp;

      // Find minimum timestep from all hydro particles
      if (hydro != NULL) {
#pragma omp for
        for (i=0; i<hydro->Nhydro; i++) {
          Particle<ndim>& part = hydro->GetParticlePointer(i);
          if (part.flags.is_dead()) continue;
          dt           = hydroint->Timestep(part,hydro);
          dt_min_aux   = min(dt_min_aux,dt);
          dt_hydro     = min(dt_hydro,dt);
          part.dt_next = dt;
        }
      }

      // Now compute minimum timestep due to stars/systems
      if (nbody != NULL) {
#pragma omp for
        for (i=0; i<nbody->Nnbody; i++) {
          dt         = nbody->Timestep(nbody->nbodydata[i]);
          dt_min_aux = min(dt_min_aux,dt);
          dt_nbody   = min(dt_nbody,dt);
          nbody->nbodydata[i]->dt_next = dt;
        }
      }

#ifdef _OPENMP
      const int ithread = omp_get_thread_num();
#else
      const int ithread = 0;
#endif

      timestep_temp[ithread]     = dt_min_aux;
      dt_min_hydro_temp[ithread] = dt_hydro;
      dt_min_nbody_temp[ithread] = dt_nbody;
    }

    for (int ithread=0; ithread<Nthreads; ithread++) {
      timestep     = min(timestep,timestep_temp[ithread]);
      dt_min_hydro = min(dt_min_hydro, dt_min_hydro_temp[ithread]);
      dt_min_nbody = min(dt_min_nbody, dt_min_nbody_temp[ithread]);
    }


    // For MPI, determine the global minimum timestep over all processors
#ifdef MPI_PARALLEL
    dt = timestep;
    MPI_Allreduce(&dt, &timestep, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = dt_min_hydro;
    MPI_Allreduce(&dt, &dt_min_hydro, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = dt_min_nbody;
    MPI_Allreduce(&dt, &dt_min_nbody, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

    // Calculate new block timestep levels.  For special case of Nlevels = 0
    // (i.e. constant timestep), then use timestep multipliers to set timestep
    if (Nlevels == 0) {
      level_max = 0;
      level_step = level_max + integration_step - 1;
      if (hydro) dt_min_hydro = simparams->floatparams["courant_mult"];
      if (nbody) dt_min_nbody = simparams->floatparams["nbody_mult"];
      dt_max = min(dt_min_hydro, dt_min_nbody);
      timestep = dt_max;
      level_max_hydro = level_max;
      level_max_nbody = level_max;
    }
    else {
      level_max  = Nlevels - 1;
      level_step = level_max + integration_step - 1;
      dt_max     = timestep*pow(2.0, level_max);
      level_max_hydro = min(ComputeTimestepLevel(dt_min_hydro, dt_max), level_max);
      level_max_nbody = min(ComputeTimestepLevel(dt_min_nbody, dt_max), level_max);
    }

    // Populate timestep levels with N-body particles.
    // Ensures that N-body particles occupy levels lower than all hydro particles
    if (nbody != NULL) {
      for (i=0; i<nbody->Nnbody; i++) {
        dt = nbody->nbodydata[i]->dt_next;
        level = min(ComputeTimestepLevel(dt, dt_max), level_max);
        nbody->nbodydata[i]->level = max(level, level_max_hydro);
        nbody->nbodydata[i]->nlast = n;
        nbody->nbodydata[i]->nstep = pow(2, level_step - nbody->nbodydata[i]->level);
        nbody->nbodydata[i]->tlast = t;
        nbody->nbodydata[i]->dt_next = nbody->nbodydata[i]->nstep * timestep;
        nbody->nbodydata[i]->flags.set(end_timestep);
      }
    }

    // Populate the timestep levels with hydro particles.
    // If particles are sink neighbours, set to same timesteps as sinks
    if (hydro != NULL) {
#pragma omp parallel default(none) private(i,dt,level) shared(level_max_nbody,level_max_temp)
      {

        int level_max_hydro_thread = 0;

#pragma omp for
        for (i=0; i<hydro->Nhydro; i++) {
          Particle<ndim>& part = hydro->GetParticlePointer(i);

          if (part.flags.is_dead()) continue;
          dt             = part.dt_next;
          level          = min(ComputeTimestepLevel(dt, dt_max), level_max);
          part.level     = level;
          part.levelneib = level;
          part.nstep     = pow(2, level_step - part.level);
          part.nlast     = n ;
          part.dt_next   = part.nstep * timestep;
          part.flags.set(end_timestep);

          if (part.flags.check(inside_sink)) {
            if (level_max_nbody - part.level > level_diff_max) {
              part.level     = level_max_nbody - level_diff_max;
              part.levelneib = level_max_nbody;
              level_max_hydro_thread = max(level_max_hydro_thread, part.level);
            }
          }
        }

#ifdef _OPENMP
        const int ithread = omp_get_thread_num();
#else
        const int ithread = 0;
#endif
        level_max_temp[ithread] = level_max_hydro_thread;
      }

      for (int ithread=0; ithread<Nthreads; ithread++) {
        level_max_hydro = max(level_max_hydro, level_max_temp[ithread]);
      }


      // If enforcing a single hydro timestep, set it here.
      // Otherwise, populate the timestep levels with hydro particles.
      if (hydro_single_timestep == 1) {
        for (i=0; i<hydro->Nhydro; i++) {
          Particle<ndim>& part = hydro->GetParticlePointer(i);

          if (part.flags.is_dead()) continue;
          dt             = part.dt_next;
          level          = min(ComputeTimestepLevel(dt, dt_max), level_max);
          part.level     = level;
          part.levelneib = level;
          part.nstep     = pow(2, level_step - part.level);
          part.nlast     = n;
          part.dt_next   = part.nstep * timestep;
          part.flags.set(end_timestep);
        }
      }
    }

    nresync = pow(2, level_step);
    assert(nresync > 0);
    timestep = dt_max / (DOUBLE) nresync;

  }
  // If not resynchronising, check if any hydro/N-body particles need to move
  // up or down timestep levels.
  //===============================================================================================
  else if (Nlevels > 0) {

    level_max_old   = level_max;
    level_max       = 0;
    level_max_nbody = 0;
    level_max_hydro = 0;


    if (hydro != NULL) {

#pragma omp parallel default(shared) private(dt,dt_nbody,dt_hydro,i)\
     private(istep,last_level,level,level_max_aux,level_nbody,level_hydro,nstep,nfactor)
      {
        dt_hydro      = big_number_dp;
        level_max_aux = 0;
        level_hydro   = 0;


        // Find all hydro particles at the beginning of a new timestep
        //-----------------------------------------------------------------------------------------
#pragma omp for
        for (i=0; i<hydro->Nhydro; i++) {
          Particle<ndim>& part = hydro->GetParticlePointer(i);
          if (part.flags.is_dead()) continue;

          // hydro particles whose timestep has been artificially reduced by Saitoh & Makino scheme.
          if (n - part.nlast == part.nstep && part.nstep != pow(2,level_step - part.level)) {
            dt             = hydroint->Timestep(part, hydro);
            level          = max(ComputeTimestepLevel(dt, dt_max), part.levelneib - level_diff_max);
            part.level     = max(part.level, level);
            part.levelneib = part.level;
            part.nlast     = n;
            part.nstep     = pow(2, level_step - part.level);
            part.dt_next   = part.nstep * timestep;
            part.flags.set(end_timestep) ;
          }
          // hydro particles that have naturally reached the end of their step
          else if (n - part.nlast == part.nstep) {
            nstep      = part.nstep;
            last_level = part.level;

            // Compute new timestep value and level number
            dt    = hydroint->Timestep(part, hydro);
            level = max(ComputeTimestepLevel(dt, dt_max), part.levelneib - level_diff_max);

            // Move up one level (if levels are correctly synchronised) or
            // down several levels if required
            if (level < last_level && last_level > 1 && n%(2*nstep) == 0) {
              part.level = last_level - 1;
            }
            else if (level > last_level) {
              part.level = level;
            }
            else {
              part.level = last_level;
            }

            part.levelneib = level;
            part.nlast     = n;
            part.nstep     = pow(2, level_step - part.level);
            part.dt_next   = part.nstep * timestep;
            part.flags.set(end_timestep);
          }

          // Find maximum level of all hydro particles
          level_hydro     = max(level_hydro, part.level);
          level_max_aux = max(level_max_aux, part.level);

          dt_hydro = min(dt_hydro, (DOUBLE) part.dt);
        }
        //-------------------------------------------------------------------------------------------

#ifdef _OPENMP
        const int ithread = omp_get_thread_num();
#else
        const int ithread = 0;
#endif
        dt_min_hydro_temp[ithread] = dt_hydro;
        level_max_temp[ithread] = level_max_aux;
      }

      for (int ithread=0; ithread<Nthreads; ithread++) {
        dt_min = min(dt_min,dt_min_hydro_temp[ithread]);
        level_max = max(level_max,level_max_temp[ithread]);
      }
      dt_min_hydro = dt_min;
      level_max_hydro = level_max;


#if defined MPI_PARALLEL
      level = level_max_hydro;
      MPI_Allreduce(&level, &level_max_hydro, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
    }

    // Now find all N-body particles at the beginning of a new timestep
    //-------------------------------------------------------------------------------------------
    if (nbody != NULL) {

      dt_nbody      = big_number_dp;
      level_max_aux = 0;
      level_nbody   = 0;

      for (i=0; i<nbody->Nnbody; i++) {

        // Skip particles that are not at end of step
        if (n - nbody->nbodydata[i]->nlast == nbody->nbodydata[i]->nstep) {
          nstep = nbody->nbodydata[i]->nstep;
          last_level = nbody->nbodydata[i]->level;

          // Compute new timestep value and level number
          dt    = nbody->Timestep(nbody->nbodydata[i]);
          level = max(ComputeTimestepLevel(dt, dt_max), level_max_hydro);

          // Move up one level (if levels are correctly synchronised) or
          // down several levels if required
          if (level < last_level && level > level_max_hydro && last_level > 1 && n%(2*nstep) == 0) {
            nbody->nbodydata[i]->level = last_level - 1;
          }
          else if (level > last_level) {
            nbody->nbodydata[i]->level = level;
          }
          else {
            nbody->nbodydata[i]->level = last_level;
          }

          nbody->nbodydata[i]->nlast = n;
          nbody->nbodydata[i]->nstep = pow(2, level_step - nbody->nbodydata[i]->level);
          nbody->nbodydata[i]->tlast = t;
          nbody->nbodydata[i]->dt_next = nbody->nbodydata[i]->nstep * timestep;
          nbody->nbodydata[i]->flags.set(end_timestep);
        }

        // Find maximum level of all N-body particles
        level_nbody   = max(level_nbody, nbody->nbodydata[i]->level);
        level_max_aux = max(level_max_aux, nbody->nbodydata[i]->level);
        dt_nbody      = min(dt_nbody, nbody->nbodydata[i]->dt_next);
      }
      //-----------------------------------------------------------------------------------------

      dt_min          = min(dt_min, dt_nbody);
      dt_min_nbody    = min(dt_min_nbody, dt_nbody);
      level_max       = max(level_max, level_max_aux);
      level_max_nbody = max(level_max_nbody, level_nbody);
    }

    // For MPI, find the global maximum timestep levels for each processor
#ifdef MPI_PARALLEL
    level = level_max;
    MPI_Allreduce(&level, &level_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    level = level_max_nbody;
    MPI_Allreduce(&level, &level_max_nbody, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    assert(level_max_hydro >= 0);
#endif

    assert(!(isnan(dt_min)) && !(isinf(dt_min)));
    assert(!(isnan(dt_max)) && !(isinf(dt_max)));


    if (hydro != NULL) {
      // Set fixed hydro timestep level here in case maximum has changed
      if (hydro_single_timestep == 1) {
        for (i=0; i<hydro->Nhydro; i++) {
          Particle<ndim>& part = hydro->GetParticlePointer(i);
          if (part.flags.is_dead()) continue;
          if (part.nlast == n) part.level = level_max_hydro;
        }
      }

      // Update all timestep variables if we have removed or added any levels
      istep = pow(2, level_step - level_max_old + 1);

      // Adjust integer time if levels are added or removed
      if (level_max > level_max_old) {
        nfactor = pow(2, level_max - level_max_old);
        n *= nfactor;

        level_step = level_max + integration_step - 1;
#pragma omp parallel for default(none) private(i) shared(nfactor)
        for (i=0; i<hydro->Nhydro; i++) {
          Particle<ndim>& part = hydro->GetParticlePointer(i);
          if (part.flags.is_dead()) continue;
          part.nstep *= nfactor;
          part.nlast *= nfactor;
          if (part.nlast == n) part.nstep = pow(2, level_step - part.level);
        }
        for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep *= nfactor;
        for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast *= nfactor;
      }
      else if (level_max <= level_max_old - 1 && level_max_old > 1 && n%istep == 0) {
        level_max = level_max_old - 1;

        level_step = level_max + integration_step - 1;

        nfactor = pow(2, level_max_old - level_max);
        assert(n%nfactor == 0);
        n /= nfactor;
#pragma omp parallel for default(none) private(i) shared(nfactor)
        for (i=0; i<hydro->Nhydro; i++) {
          Particle<ndim>& part = hydro->GetParticlePointer(i);
          if (part.flags.is_dead()) continue;
          part.nlast /= nfactor;
          part.nstep /= nfactor;
          if (part.nlast == n) part.nstep = pow(2, level_step - part.level);
        }
        for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast /= nfactor;
        for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep /= nfactor;
      }
      else {
        level_max = level_max_old;
      }
    }

    nresync    = pow(2, level_step);
    timestep   = dt_max / (DOUBLE) nresync;

    // Update values of nstep for both star particles
    if (nbody != NULL) {

      if (hydro == NULL) {
        level_step = level_max + integration_step - 1;
        nresync    = pow(2, level_step);
        timestep   = dt_max / (DOUBLE) nresync;
      }

      for (i=0; i<nbody->Nnbody; i++) {
        if (nbody->nbodydata[i]->nlast == n) {
          nbody->nbodydata[i]->nstep = pow(2, level_step - nbody->nbodydata[i]->level);
        }
      }
    }

    assert(level_max >= level_max_old - 1);

  }
  //===============================================================================================


#ifndef NDEBUG
  // Various asserts for debugging
  if (hydro != NULL) {
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      if (part.flags.is_dead()) continue;
      assert(part.level <= level_max);
      assert(part.nlast <= n);
      assert(part.nstep == pow(2,level_step - part.level));
      assert(part.nlast != n || n%part.nstep == 0);
    }
  }
  if (nbody != NULL) {
    for (i=0; i<nbody->Nnbody; i++) {
      assert(nbody->nbodydata[i]->level <= level_max);
      assert(nbody->nbodydata[i]->nlast <= n);
      assert(nbody->nbodydata[i]->nstep == pow(2,level_step - nbody->nbodydata[i]->level));
      assert(nbody->nbodydata[i]->nlast != n || n%nbody->nbodydata[i]->nstep == 0);
      assert(nbody->nbodydata[i]->level >= level_max_hydro);
      assert(nbody->nbodydata[i]->tlast <= t);
    }
  }
  assert(timestep >= 0.0 && !(isinf(timestep)) && !(isnan(timestep)));
  assert(dt_max > 0.0 && !(isinf(dt_max)) && !(isnan(dt_max)));
  assert(level_step == level_max + integration_step - 1);
  assert(level_max_hydro <= level_max);
  assert(level_max_nbody <= level_max);
  assert(n <= nresync);
  if (timestep <= 0.0) {
    cout << "Timestep fallen to zero : " << timestep << "    dtmax: " << dt_max
         << "    nresync " << nresync << endl;
    ExceptionHandler::getIstance().raise("Error : timestep fallen to zero");
  }
#endif

  return;

}



template class Simulation<1>;
template class Simulation<2>;
template class Simulation<3>;
