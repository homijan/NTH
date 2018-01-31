// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
// reserved. See files LICENSE and NOTICE for details.
//
// This file is part of CEED, a collection of benchmarks, miniapps, software
// libraries and APIs for efficient high-order finite element and spectral
// element discretizations for exascale applications. For more information and
// source code availability see http://github.com/ceed.
//
// The CEED research is supported by the Exascale Computing Project 17-SC-20-SC,
// a collaborative effort of two U.S. Department of Energy organizations (Office
// of Science and the National Nuclear Security Administration) responsible for
// the planning and preparation of a capable exascale ecosystem, including
// software, applications, hardware, advanced system engineering and early
// testbed platforms, in support of the nation's exascale computing imperative.
//
//                     __                __
//                    / /   ____  ____  / /_  ____  _____
//                   / /   / __ `/ __ `/ __ \/ __ \/ ___/
//                  / /___/ /_/ / /_/ / / / / /_/ (__  )
//                 /_____/\__,_/\__, /_/ /_/\____/____/
//                             /____/
//
//             High-order Lagrangian Hydrodynamics Miniapp
//
// Laghos(LAGrangian High-Order Solver) is a miniapp that solves the
// time-dependent Euler equation of compressible gas dynamics in a moving
// Lagrangian frame using unstructured high-order finite element spatial
// discretization and explicit high-order time-stepping. Laghos is based on the
// numerical algorithm described in the following article:
//
//    V. Dobrev, Tz. Kolev and R. Rieben, "High-order curvilinear finite element
//    methods for Lagrangian hydrodynamics", SIAM Journal on Scientific
//    Computing, (34) 2012, pp.B606–B641, https://doi.org/10.1137/120864672.
//
// Sample runs:
//    mpirun -np 8 laghos -p 0 -m data/square01_quad.mesh -rs 3 -tf 0.75
//    mpirun -np 8 laghos -p 0 -m data/square01_tri.mesh  -rs 1 -tf 0.75
//    mpirun -np 8 laghos -p 0 -m data/cube01_hex.mesh    -rs 1 -tf 2.0
//    mpirun -np 8 laghos -p 1 -m data/square01_quad.mesh -rs 3 -tf 0.8
//    mpirun -np 8 laghos -p 1 -m data/square01_quad.mesh -rs 0 -tf 0.8 -ok 7 -ot 6
//    mpirun -np 8 laghos -p 1 -m data/cube01_hex.mesh    -rs 2 -tf 0.6
//    mpirun -np 8 laghos -p 2 -m data/segment01.mesh     -rs 5 -tf 0.2
//    mpirun -np 8 laghos -p 3 -m data/rectangle01_quad.mesh -rs 2 -tf 2.5
//    mpirun -np 8 laghos -p 3 -m data/box01_hex.mesh        -rs 1 -tf 2.5
//
// Test problems:
//    p = 0  --> Taylor-Green vortex (smooth problem).
//    p = 1  --> Sedov blast.
//    p = 2  --> 1D Sod shock tube.
//    p = 3  --> Triple point.


#include "laghos_solver.hpp"
#include "m1_solver.hpp"
#include "eos.hpp"
#include "ic.hpp"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;
using namespace mfem::hydrodynamics;

// Choice for the problem setup.
// int problem; 

void display_banner(ostream & os);

int main(int argc, char *argv[])
{
   // Initialize MPI.
   MPI_Session mpi(argc, argv);
   int myid = mpi.WorldRank();

   // Print the banner.
   if (mpi.Root()) { display_banner(cout); }

   // Parse command-line options.
   const char *mesh_file = "data/square01_quad.mesh";
   int rs_levels = 0;
   int rp_levels = 0;
   int order_v = 2;
   int order_e = 1;
   int ode_solver_type = 4;
   double t_final = 0.5;
   double cfl = 0.5;
   double cg_tol = 1e-8;
   int cg_max_iter = 300;
   int max_tsteps = -1;
   bool p_assembly = true;
   bool visualization = false;
   int vis_steps = 10;
   bool visit = false;
   bool gfprint = false;
   const char *basename = "results/M1hos";
   int nth_problem = 5;
   double T_max = 1000.0, T_min = 100.0, rho_max = 10.0, rho_min = 1.0;
   double T_gradscale = 50.0, rho_gradscale = 50.0;
   double a0 = 1e20;
   double Zbar = 47.0;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&rs_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&rp_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&nth_problem, "-p", "--problem", "Problem setup to use.");
   args.AddOption(&order_v, "-ok", "--order-kinematic",
                  "Order (degree) of the kinematic finite element space.");
   args.AddOption(&order_e, "-ot", "--order-thermo",
                  "Order (degree) of the thermodynamic finite element space.");
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  "ODE solver: 1 - Forward Euler,\n\t"
                  "            2 - RK2 SSP, 3 - RK3 SSP, 4 - RK4, 6 - RK6.");
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&cfl, "-cfl", "--cfl", "CFL-condition number.");
   args.AddOption(&cg_tol, "-cgt", "--cg-tol",
                  "Relative CG tolerance (velocity linear solve).");
   args.AddOption(&cg_max_iter, "-cgm", "--cg-max-steps",
                  "Maximum number of CG iterations (velocity linear solve).");
   args.AddOption(&max_tsteps, "-ms", "--max-steps",
                  "Maximum number of steps (negative means no restriction).");
   args.AddOption(&p_assembly, "-pa", "--partial-assembly", "-fa",
                  "--full-assembly",
                  "Activate 1D tensor-based assembly (partial assembly).");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.AddOption(&visit, "-visit", "--visit", "-no-visit", "--no-visit",
                  "Enable or disable VisIt visualization.");
   args.AddOption(&gfprint, "-print", "--print", "-no-print", "--no-print",
                  "Enable or disable result output (files in mfem format).");
   args.AddOption(&basename, "-k", "--outputfilename",
                  "Name of the visit dump files");
   args.AddOption(&a0, "-a0", "--a0",
                  "Mean-free-path scaling, i.e. lambda = v^4/rho/a0.");
   args.AddOption(&Zbar, "-Z", "--Zbar",
                  "Constant ionization used for nu_ee. Used along IGEOS only.");
   args.AddOption(&T_max, "-Tmax", "--Tmax",
                  "Maximum temperature in the step function tanh(x).");
   args.AddOption(&T_min, "-Tmin", "--Tmin",
                  "Minimum temperature in the step function tanh(x).");
   args.AddOption(&rho_max, "-rmax", "--rhomax",
                  "Maximum density in the step function tanh(x).");
   args.AddOption(&rho_min, "-rmin", "--rhomin",
                  "Minimum density in the step function tanh(x).");
   args.AddOption(&T_gradscale, "-Tgrad", "--Tgrad",
                  "Temperature gradient scale in the function tanh(a*x).");
   args.AddOption(&rho_gradscale, "-rgrad", "--rhograd",
                  "Density gradient scale in the function tanh(a*x).");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.Root()) { args.PrintUsage(cout); }
      return 1;
   }
   if (mpi.Root()) { args.PrintOptions(cout); }

   nth::nth_problem = nth_problem;
   nth::T_max = T_max;
   nth::T_min = T_min;
   nth::rho_max = rho_max;
   nth::rho_min = rho_min;
   nth::T_gradscale = T_gradscale;
   nth::rho_gradscale = rho_gradscale;
   nth::a0 = a0;

   // Read the serial mesh from the given mesh file on all processors.
   // Refine the mesh in serial to increase the resolution.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   const int dim = mesh->Dimension();
   for (int lev = 0; lev < rs_levels; lev++) { mesh->UniformRefinement(); }

   if (p_assembly && dim == 1)
   {
      p_assembly = false;
      if (mpi.Root())
      {
         cout << "Laghos does not support PA in 1D. Switching to FA." << endl;
      }
   }

   // Parallel partitioning of the mesh.
   ParMesh *pmesh = NULL;
   const int num_tasks = mpi.WorldSize();
   const int partitions = floor(pow(num_tasks, 1.0 / dim) + 1e-2);
   int *nxyz = new int[dim];
   int product = 1;
   for (int d = 0; d < dim; d++)
   {
      nxyz[d] = partitions;
      product *= partitions;
   }
   if (product == num_tasks)
   {
      int *partitioning = mesh->CartesianPartitioning(nxyz);
      pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, partitioning);
      delete partitioning;
   }
   else
   {
      if (myid == 0)
      {
         cout << "Non-Cartesian partitioning through METIS will be used.\n";
#ifndef MFEM_USE_METIS
         cout << "MFEM was built without METIS. "
              << "Adjust the number of tasks to use a Cartesian split." << endl;
#endif
      }
#ifndef MFEM_USE_METIS
      return 1;
#endif
      pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   }
   delete [] nxyz;
   delete mesh;

   // Refine the mesh further in parallel to increase the resolution.
   for (int lev = 0; lev < rp_levels; lev++) { pmesh->UniformRefinement(); }

   // Define the parallel finite element spaces. We use:
   // - H1 (Gauss-Lobatto, continuous) for position and velocity.
   // - L2 (Bernstein, discontinuous) for specific internal energy.
   L2_FECollection L2FEC(order_e, dim, BasisType::Positive);
   H1_FECollection H1FEC(order_v, dim);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, pmesh->Dimension());

   // Boundary conditions: all tests use v.n = 0 on the boundary, and we assume
   // that the boundaries are straight.
   Array<int> ess_tdofs;
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max()), tdofs1d;
      for (int d = 0; d < pmesh->Dimension(); d++)
      {
         // Attributes 1/2/3 correspond to fixed-x/y/z boundaries, i.e., we must
         // enforce v_x/y/z = 0 for the velocity components.
         ess_bdr = 0; ess_bdr[d] = 1;
         H1FESpace.GetEssentialTrueDofs(ess_bdr, tdofs1d, d);
         ess_tdofs.Append(tdofs1d);
      }
   }

   // Define the explicit ODE solver used for time integration.
   ODESolver *ode_solver = NULL;
   switch (ode_solver_type)
   {
      case 1: ode_solver = new ForwardEulerSolver; break;
      case 2: ode_solver = new RK2Solver(0.5); break;
      case 3: ode_solver = new RK3SSPSolver; break;
      case 4: ode_solver = new RK4Solver; break;
      case 6: ode_solver = new RK6Solver; break;
      default:
         if (myid == 0)
         {
            cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
         }
         delete pmesh;
         MPI_Finalize();
         return 3;
   }

   HYPRE_Int glob_size_l2 = L2FESpace.GlobalTrueVSize();
   HYPRE_Int glob_size_h1 = H1FESpace.GlobalTrueVSize();

   if (mpi.Root())
   {
      cout << "Number of kinematic (position, velocity) dofs: "
           << glob_size_h1 << endl;
      cout << "Number of specific internal energy dofs: "
           << glob_size_l2 << endl;
   }

   int Vsize_l2 = L2FESpace.GetVSize();
   int Vsize_h1 = H1FESpace.GetVSize();

   // The monolithic BlockVector stores unknown fields as:
   // - 0 -> position
   // - 1 -> velocity
   // - 2 -> specific internal energy

   Array<int> true_offset(4);
   true_offset[0] = 0;
   true_offset[1] = true_offset[0] + Vsize_h1;
   true_offset[2] = true_offset[1] + Vsize_h1;
   true_offset[3] = true_offset[2] + Vsize_l2;
   BlockVector S(true_offset);

   // Define GridFunction objects for the position, velocity and specific
   // internal energy.  There is no function for the density, as we can always
   // compute the density values given the current mesh position, using the
   // property of pointwise mass conservation.
   ParGridFunction x_gf, v_gf, e_gf;
   x_gf.MakeRef(&H1FESpace, S, true_offset[0]);
   v_gf.MakeRef(&H1FESpace, S, true_offset[1]);
   e_gf.MakeRef(&L2FESpace, S, true_offset[2]);

   // Initialize x_gf using the starting mesh coordinates. This also links the
   // mesh positions to the values in x_gf.
   pmesh->SetNodalGridFunction(&x_gf);

   // Initialize the velocity.
   VectorFunctionCoefficient v_coeff(pmesh->Dimension(), nth::v0);
   v_gf.ProjectCoefficient(v_coeff);

   // Initialize density and specific internal energy values. We interpolate in
   // a non-positive basis to get the correct values at the dofs.  Then we do an
   // L2 projection to the positive basis in which we actually compute. The goal
   // is to get a high-order representation of the initial condition. Note that
   // this density is a temporary function and it will not be updated during the
   // time evolution.
   ParGridFunction rho_gf(&L2FESpace);
   FunctionCoefficient rho_coeff(nth::rho0);
   L2_FECollection l2_fec(order_e, pmesh->Dimension());
   ParFiniteElementSpace l2_fes(pmesh, &l2_fec);
   ParGridFunction l2_rho(&l2_fes), l2_e(&l2_fes);
   l2_rho.ProjectCoefficient(rho_coeff);
   rho_gf.ProjectGridFunction(l2_rho);
   if (nth::nth_problem == 1)
   {
      // For the Sedov test, we use a delta function at the origin.
      DeltaCoefficient e_coeff(0, 0, 0.25);
      l2_e.ProjectCoefficient(e_coeff);
	  // Set min temperature to be nonzero.
      //ParGridFunction cnst(&l2_fes);
	  //cnst = 0.0025;
	  //l2_e += cnst;
   }
   else
   {
      FunctionCoefficient e_coeff(nth::e0);
      l2_e.ProjectCoefficient(e_coeff);
   }
   e_gf.ProjectGridFunction(l2_e);

   // Space-dependent ideal gas coefficient over the Lagrangian mesh.
   FunctionCoefficient gamma_cf = (nth::gamma);
   Coefficient *material_pcf = &gamma_cf;

   // Additional details, depending on the problem.
   int source = 0; bool visc;
   switch (nth::nth_problem)
   {
      case 0: if (pmesh->Dimension() == 2) { source = 1; }
         visc = false; break;
      case 1: visc = true; break;
      case 2: visc = true; break;
      case 3: visc = true; break;
      case 4: visc = true; break;
      case 5: visc = true; break;
      case 6: visc = true; break;
      case 7: visc = true; break;
      case 8: visc = true; break;
      default: MFEM_ABORT("Wrong problem specification!");
   }

   LagrangianHydroOperator oper(S.Size(), H1FESpace, L2FESpace,
                                ess_tdofs, rho_gf, source, cfl, material_pcf,
                                visc, p_assembly, cg_tol, cg_max_iter);

///////////////////////////////////////////////////////////////
///// M1 nonlocal solver //////////////////////////////////////
///////////////////////////////////////////////////////////////
   // The monolithic BlockVector stores unknown fields as:
   // - 0 -> isotropic I0 (energy density)
   // - 1 -> anisotropic I1 (flux density)
   Array<int> m1true_offset(3);
   m1true_offset[0] = 0;
   m1true_offset[1] = m1true_offset[0] + Vsize_l2;
   m1true_offset[2] = m1true_offset[1] + Vsize_h1;
   BlockVector m1S(m1true_offset);

   // Define GridFunction objects for the zero and first moments of
   // the~electron distribution function.
   ParGridFunction I0_gf, I1_gf;
   I0_gf.MakeRef(&L2FESpace, m1S, m1true_offset[0]);
   I1_gf.MakeRef(&H1FESpace, m1S, m1true_offset[1]);

   // Define hydrodynamics related coefficients as mean stopping power and
   // source function depending on plasma temperature and density. 
   const double kB = 1.0, me = 1.0, pi = 3.14159265359;
   nth::IGEOS eos(me, kB);
   // Use a constant ionization provided by IG eos. 
   eos.SetZbar(Zbar);
   // Prepare C6 physics.
   nth::ClassicalMeanStoppingPower mspei_cf(rho_gf, e_gf, v_gf, material_pcf,
                                            &eos);
   nth::ClassicalAWBSMeanStoppingPower mspee_cf(rho_gf, e_gf, v_gf, 
                                                material_pcf, &eos);
   nth::NTHvHydroCoefficient *mspei_pcf = &mspei_cf;
   nth::NTHvHydroCoefficient *mspee_pcf = &mspee_cf;
   nth::ClassicalMeanFreePath mfp_cf(rho_gf, e_gf, v_gf, material_pcf, &eos);
   nth::MeanFreePath *mfp_pcf = &mfp_cf;
   nth::AWBSI0Source sourceI0_cf(rho_gf, e_gf, v_gf, material_pcf, &eos);
   nth::NTHvHydroCoefficient *sourceI0_pcf = &sourceI0_cf;
   nth::KnudsenNumber Kn_cf(rho_gf, e_gf, v_gf, material_pcf, &eos, mfp_pcf);
   Coefficient *Kn_pcf = &Kn_cf;
   nth::LorentzEfield LorentzEfield_cf(pmesh->Dimension(), rho_gf, e_gf, v_gf, 
                                material_pcf, &eos);
   VectorCoefficient *Efield_pcf = &LorentzEfield_cf;
   Vector vZero(pmesh->Dimension());
   vZero = 0.0;
   VectorConstantCoefficient ZeroBfield_cf(vZero);
   VectorCoefficient *Bfield_pcf = &ZeroBfield_cf;

   nth::AWBSMasterOfPhysics AWBSPhysics(mspei_pcf, mspee_pcf, sourceI0_pcf,
                                        Efield_pcf, Bfield_pcf);

   // Static coefficient defined in m1_solver.hpp.
   double m1cfl = 0.25;
   vis_steps = 1000000000;
   // ALWAYS calculate on v in (0, 1)
   double vmax = 1.0;
   double vTmultiple = 7.0;
   // well, not really, since the lowest v = 0 is singular, so
   //double vmin = 0.01 * vmax;
   double vmin = 0.02 * vmax;
   // and provide some maximum dv step.
   double dvmax = vmax*0.0005;
   bool nonlocal_test = false;
   if (nonlocal_test)
   {
      vTmultiple = 6.0;
      vmin = 0.01 * vmax;
      //vmin = 3.5 * vmax / vTmultiple; // Minimum 3.5*vTh
      dvmax = vmax*0.1;
      if (pmesh->Dimension() == 1)
      { 
         nth::a0 = 2e3;
         vis_steps = 10000;
         m1cfl = 0.5;
      }
      else if (pmesh->Dimension() == 2)
      { 
	     nth::a0 = 1e5; //2e1;
         vis_steps = 10000;
         m1cfl = 1.0;
      }
      else if (pmesh->Dimension() == 3)
      {
         nth::a0 = 5e7;
         vis_steps = 10000;
         m1cfl = 0.5;
      }
   }

   oper.ComputeDensity(rho_gf);
   AWBSPhysics.SetThermalVelocityMultiple(vTmultiple);
   mfp_cf.SetThermalVelocityMultiple(vTmultiple);
   double loc_Tmax = e_gf.Max(), glob_Tmax;
   MPI_Allreduce(&loc_Tmax, &glob_Tmax, 1, MPI_DOUBLE, MPI_MAX,
                 pmesh->GetComm());
   AWBSPhysics.SetTmax(glob_Tmax);
   mfp_cf.SetTmax(glob_Tmax);

   // Initialize the M1-AWBS operator
   nth::M1Operator m1oper(m1S.Size(), H1FESpace, L2FESpace, ess_tdofs, rho_gf, 
                          m1cfl, &AWBSPhysics, x_gf, e_gf, cg_tol, cg_max_iter);
   // Prepare grid functions integrating the moments of I0 and I1.
   ParGridFunction intf0_gf(&L2FESpace), Kn_gf(&L2FESpace);
   ParGridFunction j_gf(&H1FESpace), hflux_gf(&H1FESpace);

   ODESolver *m1ode_solver = NULL;
   //m1ode_solver = new ForwardEulerSolver;
   //m1ode_solver = new RK2Solver(0.5);
   m1ode_solver = new RK4Solver;
   //m1ode_solver = new RK6Solver;
   m1ode_solver->Init(m1oper);

   double alphavT = mspei_cf.GetVelocityScale();
   m1oper.ResetVelocityStepEstimate();
   m1oper.ResetQuadratureData();
   m1oper.SetTime(vmax);
   //double dvmax = vmax*0.1;
   double dvmin = min(dvmax, m1oper.GetVelocityStepEstimate(m1S));
   I0_gf = 0.0; I1_gf = 0.0;
   int m1ti = 0;
   double v = vmax;
   double dv = -dvmin;
   intf0_gf = 0.0;
   j_gf = 0.0;
   hflux_gf = 0.0;
   Kn_gf.ProjectCoefficient(Kn_cf);
/*
   while (abs(dv) >= abs(dvmin))
   {
      m1ti++;
      m1ode_solver->Step(m1S, v, dv);

      // Perform the integration over velocity space.
      intf0_gf.Add(pow(alphavT*v, 2.0) * alphavT*abs(dv), I0_gf);
      j_gf.Add(pow(alphavT*v, 3.0) * alphavT*abs(dv), I1_gf);
      hflux_gf.Add(me / 2.0 * pow(alphavT*v, 5.0) * alphavT*abs(dv), I1_gf);

      double loc_minI0 = I0_gf.Min(), glob_minI0;
      MPI_Allreduce(&loc_minI0, &glob_minI0, 1, MPI_DOUBLE, MPI_MIN,
                       pmesh->GetComm());
      double loc_maxI0 = I0_gf.Max(), glob_maxI0;
      MPI_Allreduce(&loc_maxI0, &glob_maxI0, 1, MPI_DOUBLE, MPI_MAX,
                       pmesh->GetComm());
      double loc_minI1 = I1_gf.Min(), glob_minI1;
      MPI_Allreduce(&loc_minI1, &glob_minI1, 1, MPI_DOUBLE, MPI_MIN,
                       pmesh->GetComm());
      double loc_maxI1 = I1_gf.Max(), glob_maxI1;
      MPI_Allreduce(&loc_maxI1, &glob_maxI1, 1, MPI_DOUBLE, MPI_MAX,
                       pmesh->GetComm());

      m1oper.ResetVelocityStepEstimate();
      m1oper.ResetQuadratureData();
      m1oper.SetTime(v);
      dv = - min(dvmax, m1oper.GetVelocityStepEstimate(m1S));
      if (v + dv < vmin) { dv = vmin - v; }

      if (mpi.Root())
      {
         cout << fixed;
         cout << "group " << setw(5) << m1ti
                 << ",\tv = " << setw(5) << setprecision(4) << v
                 << ",\tdv = " << setw(5) << setprecision(8) << dv << endl
                 << "[min(f0), max(f0)] = [" << setprecision(17)
                 << glob_minI0 << ",\t" << glob_maxI0 << "]" << endl
                 << "[min(f1), max(f1)] = [" << setprecision(17)
                 << glob_minI1 << ",\t" << glob_maxI1 << "]"
                 << endl;
      }
   }
*/
///////////////////////////////////////////////////////////////
///// M1 nonlocal solver //////////////////////////////////////
///////////////////////////////////////////////////////////////

   socketstream vis_rho, vis_v, vis_e, vis_f0, vis_j, vis_Kn, vis_hflux;
   char vishost[] = "localhost";
   int  visport   = 19916;

   if (visualization || visit) { oper.ComputeDensity(rho_gf); }

   if (visualization)
   {
      // Make sure all MPI ranks have sent their 'v' solution before initiating
      // another set of GLVis connections (one from each rank):
      MPI_Barrier(pmesh->GetComm());

      vis_rho.precision(8);
      vis_v.precision(8);
      vis_e.precision(8);

      vis_f0.precision(8);
      vis_j.precision(8);
      vis_Kn.precision(8);
      vis_hflux.precision(8);

      int Wx = 0, Wy = 0; // window position
      const int Ww = 350, Wh = 350; // window size
      int offx = Ww+10; // window offsets

      VisualizeField(vis_rho, vishost, visport, rho_gf,
                     "Density", Wx, Wy, Ww, Wh);
      //Wx += offx;
      //VisualizeField(vis_v, vishost, visport, v_gf,
      //               "Velocity", Wx, Wy, Ww, Wh);
      Wx += offx;
      VisualizeField(vis_j, vishost, visport, j_gf,
                     "Current", Wx, Wy, Ww, Wh);
      Wx += offx;
      VisualizeField(vis_e, vishost, visport, e_gf,
                     "T", Wx, Wy, Ww, Wh);

      Wx = 0;
      Wy +=offx;
      VisualizeField(vis_f0, vishost, visport, intf0_gf,
                     "int(f0 4pi v^2)dv", Wx, Wy, Ww, Wh);
      Wx += offx;
      VisualizeField(vis_Kn, vishost, visport, Kn_gf,
                     "Kn", Wx, Wy, Ww, Wh);
      Wx += offx;
      VisualizeField(vis_hflux, vishost, visport, hflux_gf,
                     "Heat flux", Wx, Wy, Ww, Wh);

   }

   // Save data for VisIt visualization
   VisItDataCollection visit_dc(basename, pmesh);
   if (visit)
   {
      visit_dc.RegisterField("Density",  &rho_gf);
      visit_dc.RegisterField("Velocity", &v_gf);
      visit_dc.RegisterField("Specific Internal Energy", &e_gf);
      visit_dc.SetCycle(0);
      visit_dc.SetTime(0.0);
      visit_dc.Save();
   }

   // Perform time-integration (looping over the time iterations, ti, with a
   // time-step dt). The object oper is of type LagrangianHydroOperator that
   // defines the Mult() method that used by the time integrators.
   ode_solver->Init(oper);
   oper.ResetTimeStepEstimate();
   double t = 0.0, dt = oper.GetTimeStepEstimate(S), t_old;
   bool last_step = false;
   int steps = 0;
   BlockVector S_old(S);
   for (int ti = 1; !last_step; ti++)
   {
      if (t + dt >= t_final)
      {
         dt = t_final - t;
         last_step = true;
      }
      if (steps == max_tsteps) { last_step = true; }

      S_old = S;
      t_old = t;
      oper.ResetTimeStepEstimate();

      // S is the vector of dofs, t is the current time, and dt is the time step
      // to advance.
      ode_solver->Step(S, t, dt);
      steps++;

      // Adaptive time step control.
      const double dt_est = oper.GetTimeStepEstimate(S);
      if (dt_est < dt)
      {
         // Repeat (solve again) with a decreased time step - decrease of the
         // time estimate suggests appearance of oscillations.
         dt *= 0.85;
         if (dt < numeric_limits<double>::epsilon())
         { MFEM_ABORT("The time step crashed!"); }
         t = t_old;
         S = S_old;
         oper.ResetQuadratureData();
         if (mpi.Root()) { cout << "Repeating step " << ti << endl; }
         ti--; continue;
      }
      else if (dt_est > 1.25 * dt) { dt *= 1.02; }

      // Make sure that the mesh corresponds to the new solution state.
      pmesh->NewNodes(x_gf, false);

      if (last_step || (ti % vis_steps) == 0)
      {
         double loc_norm = e_gf * e_gf, tot_norm;
         MPI_Allreduce(&loc_norm, &tot_norm, 1, MPI_DOUBLE, MPI_SUM,
                       pmesh->GetComm());
         if (mpi.Root())
         {
            cout << fixed;
            cout << "step " << setw(5) << ti
                 << ",\tt = " << setw(5) << setprecision(4) << t
                 << ",\tdt = " << setw(5) << setprecision(6) << dt
                 << ",\t|e| = " << setprecision(10)
                 << sqrt(tot_norm) << endl;
         }

         // Make sure all ranks have sent their 'v' solution before initiating
         // another set of GLVis connections (one from each rank):
         MPI_Barrier(pmesh->GetComm());

///////////////////////////////////////////////////////////////
///// M1 nonlocal solver //////////////////////////////////////
///////////////////////////////////////////////////////////////
         oper.ComputeDensity(rho_gf);
         AWBSPhysics.SetThermalVelocityMultiple(vTmultiple);
		 mfp_cf.SetThermalVelocityMultiple(vTmultiple);          
         double loc_Tmax = e_gf.Max(), glob_Tmax;
         MPI_Allreduce(&loc_Tmax, &glob_Tmax, 1, MPI_DOUBLE, MPI_MAX,
                       pmesh->GetComm());
         AWBSPhysics.SetTmax(glob_Tmax);
		 mfp_cf.SetTmax(glob_Tmax);
         alphavT = mspei_cf.GetVelocityScale();
         m1oper.ResetVelocityStepEstimate();
         m1oper.ResetQuadratureData();
         m1oper.SetTime(vmax);
         double dvmin = min(dvmax, m1oper.GetVelocityStepEstimate(m1S));
         I0_gf = 0.0; //1e-2; 
		 I1_gf = 0.0;
         int m1ti = 0;
         double v = vmax;
         double dv = -dvmin;
         intf0_gf = 0.0;
         j_gf = 0.0;
         hflux_gf = 0.0;
         Kn_gf.ProjectCoefficient(Kn_cf);
         while (abs(dv) >= abs(dvmin))
         {
            m1ti++;
            m1ode_solver->Step(m1S, v, dv);

            // Perform the integration over velocity space.
            intf0_gf.Add(pow(alphavT*v, 2.0) * alphavT*abs(dv), I0_gf);
            j_gf.Add(pow(alphavT*v, 3.0) * alphavT*abs(dv), I1_gf);
            hflux_gf.Add(me / 2.0 * pow(alphavT*v, 5.0) * alphavT*abs(dv), 
                         I1_gf);

			double loc_minI0 = I0_gf.Min(), glob_minI0;
            MPI_Allreduce(&loc_minI0, &glob_minI0, 1, MPI_DOUBLE, MPI_MIN,
                          pmesh->GetComm());
            double loc_maxI0 = I0_gf.Max(), glob_maxI0;
            MPI_Allreduce(&loc_maxI0, &glob_maxI0, 1, MPI_DOUBLE, MPI_MAX,
                          pmesh->GetComm());
            double loc_minI1 = I1_gf.Min(), glob_minI1;
            MPI_Allreduce(&loc_minI1, &glob_minI1, 1, MPI_DOUBLE, MPI_MIN,
                          pmesh->GetComm());
            double loc_maxI1 = I1_gf.Max(), glob_maxI1;
            MPI_Allreduce(&loc_maxI1, &glob_maxI1, 1, MPI_DOUBLE, MPI_MAX,
                          pmesh->GetComm());

            m1oper.ResetVelocityStepEstimate();
            m1oper.ResetQuadratureData();
            m1oper.SetTime(v);
            dv = - min(dvmax, m1oper.GetVelocityStepEstimate(m1S));
            if (v + dv < vmin) { dv = vmin - v; }

            if (mpi.Root())
            {
               cout << fixed;
               cout << "group " << setw(5) << m1ti
               << ",\tv = " << setw(5) << setprecision(4) << v
               << ",\tdv = " << setw(5) << setprecision(8) << dv << endl
               << "[min(f0), max(f0)] = [" << setprecision(17)
               << glob_minI0 << ",\t" << glob_maxI0 << "]" << endl
               << "[min(f1), max(f1)] = [" << setprecision(17)
               << glob_minI1 << ",\t" << glob_maxI1 << "]"
               << endl;
            }
         }
///////////////////////////////////////////////////////////////
///// M1 nonlocal solver //////////////////////////////////////
///////////////////////////////////////////////////////////////

         if (visualization || visit || gfprint) { oper.ComputeDensity(rho_gf); }
         if (visualization)
         {
            int Wx = 0, Wy = 0; // window position
            int Ww = 350, Wh = 350; // window size
            int offx = Ww+10; // window offsets

            VisualizeField(vis_rho, vishost, visport, rho_gf,
                           "Density", Wx, Wy, Ww, Wh);
            //Wx += offx;
            //VisualizeField(vis_v, vishost, visport,
            //               v_gf, "Velocity", Wx, Wy, Ww, Wh);
            Wx += offx;
            VisualizeField(vis_j, vishost, visport, j_gf,
                           "Current", Wx, Wy, Ww, Wh);		
            Wx += offx;
            VisualizeField(vis_e, vishost, visport, e_gf,
                           "T", Wx, Wy, Ww,Wh);

            Wx = 0;
            Wy +=offx;
            VisualizeField(vis_f0, vishost, visport, intf0_gf,
                           "int(f0 4pi v^2)dv", Wx, Wy, Ww, Wh);
            Wx += offx;
            VisualizeField(vis_Kn, vishost, visport, Kn_gf,
                           "Kn", Wx, Wy, Ww, Wh);
            Wx += offx;
            VisualizeField(vis_hflux, vishost, visport, hflux_gf,
                           "Heat flux", Wx, Wy, Ww, Wh);
         }

         if (visit)
         {
            visit_dc.SetCycle(ti);
            visit_dc.SetTime(t);
            visit_dc.Save();
         }

         if (gfprint)
         {
            ostringstream mesh_name, rho_name, v_name, e_name;
            mesh_name << basename << "_" << ti
                      << "_mesh." << setfill('0') << setw(6) << myid;
            rho_name  << basename << "_" << ti
                      << "_rho." << setfill('0') << setw(6) << myid;
            v_name << basename << "_" << ti
                   << "_v." << setfill('0') << setw(6) << myid;
            e_name << basename << "_" << ti
                   << "_e." << setfill('0') << setw(6) << myid;

            ofstream mesh_ofs(mesh_name.str().c_str());
            mesh_ofs.precision(8);
            pmesh->Print(mesh_ofs);
            mesh_ofs.close();

            ofstream rho_ofs(rho_name.str().c_str());
            rho_ofs.precision(8);
            rho_gf.Save(rho_ofs);
            rho_ofs.close();

            ofstream v_ofs(v_name.str().c_str());
            v_ofs.precision(8);
            v_gf.Save(v_ofs);
            v_ofs.close();

            ofstream e_ofs(e_name.str().c_str());
            e_ofs.precision(8);
            e_gf.Save(e_ofs);
            e_ofs.close();

            ostringstream intf0_name, j_name, Kn_name, hflux_name;
            intf0_name << basename << "_" << ti
                   << "_f0." << setfill('0') << setw(6) << myid;
            j_name << basename << "_" << ti
                   << "_j." << setfill('0') << setw(6) << myid;
            Kn_name << basename << "_" << ti
                   << "_Kn." << setfill('0') << setw(6) << myid;
            hflux_name << basename << "_" << ti
                   << "_hflux." << setfill('0') << setw(6) << myid;

            ofstream intf0_ofs(intf0_name.str().c_str());
            intf0_ofs.precision(8);
            intf0_gf.Save(intf0_ofs);
            intf0_ofs.close();

            ofstream j_ofs(j_name.str().c_str());
            j_ofs.precision(8);
            j_gf.Save(j_ofs);
            j_ofs.close();

            ofstream Kn_ofs(Kn_name.str().c_str());
            Kn_ofs.precision(8);
            Kn_gf.Save(Kn_ofs);
            Kn_ofs.close();

            ofstream hflux_ofs(hflux_name.str().c_str());
            hflux_ofs.precision(8);
            hflux_gf.Save(hflux_ofs);
            hflux_ofs.close();
         }
      }
   }

   switch (ode_solver_type)
   {
      case 2: steps *= 2; break;
      case 3: steps *= 3; break;
      case 4: steps *= 4; break;
      case 6: steps *= 6;
   }
   if (mpi.Root()) { cout << "Hydrodynamics kernel timer:" << endl << flush; }
   oper.PrintTimingData(mpi.Root(), steps);
   if (mpi.Root()) { cout << "M1 kernel timer:" << endl << flush; }
   m1oper.PrintTimingData(mpi.Root(), steps);

   if (visualization)
   {
      vis_v.close();
      vis_e.close();
   }

   // Free the used memory.
   delete ode_solver;
   delete pmesh;
   delete tensors1D; 
   delete m1ode_solver;
   delete nth::tensors1D;

   return 0;
}

void display_banner(ostream & os)
{
   os << endl
      << "       __                __                 " << endl
      << "      / /   ____  ____  / /_  ____  _____   " << endl
      << "     / /   / __ `/ __ `/ __ \\/ __ \\/ ___/ " << endl
      << "    / /___/ /_/ / /_/ / / / / /_/ (__  )    " << endl
      << "   /_____/\\__,_/\\__, /_/ /_/\\____/____/  " << endl
      << "               /____/                       " << endl << endl;
}
