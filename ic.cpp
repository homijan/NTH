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

#include "ic.hpp"

#ifdef MFEM_USE_MPI

using namespace std;

namespace mfem
{

namespace nth
{

int nth_problem = 1;

} // namespace nth

namespace hydrodynamics
{

int hydro_problem = 1;
double T_max = 1000.0, T_min = 100.0;
double rho_max = 10.0, rho_min = 1.0;
double T_gradscale = 50.0, rho_gradscale = 50.0;

double rho0(const Vector &x)
{
   switch (hydro_problem)
   {
      case 0: return 1.0;
      case 1: return 1.0;
      case 2: if (x(0) < 0.5) { return 1.0; }
         else { return 0.1; }
      case 3: if (x(0) > 1.0 && x(1) <= 1.5) { return 1.0; }
         else { return 0.125; }
      case 4: return rho_min + (rho_max - rho_min) *
                     (0.5 * (tanh(rho_gradscale * (x.Norml2() - 0.5)) + 1.0));
      case 5: return 1.0;
      case 6: return rho_min + (rho_max - rho_min) *
                     (0.5 * (tanh(rho_gradscale * (x.Norml2() - 0.5)) + 1.0));
      default: MFEM_ABORT("Bad number given for problem id!"); return 0.0;
   }
}

double gamma(const Vector &x)
{
   switch (hydro_problem)
   {
      case 0: return 5./3.;
      case 1: return 1.4;
      case 2: return 1.4;
      case 3: if (x(0) > 1.0 && x(1) <= 1.5) { return 1.4; }
         else { return 1.5; }
      case 4: return 1.4;
      case 5: return 1.4;
      case 6: return 1.4;
      default: MFEM_ABORT("Bad number given for problem id!"); return 0.0;
   }
}

void v0(const Vector &x, Vector &v)
{
   switch (hydro_problem)
   {
      case 0:
         v(0) =  sin(M_PI*x(0)) * cos(M_PI*x(1));
         v(1) = -cos(M_PI*x(0)) * sin(M_PI*x(1));
         if (x.Size() == 3)
         {
            v(0) *= cos(M_PI*x(2));
            v(1) *= cos(M_PI*x(2));
            v(2) = 0.0;
         }
         break;
      case 1: v = 0.0; break;
      case 2: v = 0.0; break;
      case 3: v = 0.0; break;
      case 4: v = 0.0; break;
      case 5: v = 0.0; break;
      case 6: v = 0.0; break;
      default: MFEM_ABORT("Bad number given for problem id!");
   }
}

double e0(const Vector &x)
{
   switch (hydro_problem)
   {
      case 0:
      {
         const double denom = 2.0 / 3.0;  // (5/3 - 1) * density.
         double val;
         if (x.Size() == 2)
         {
            val = 1.0 + (cos(2*M_PI*x(0)) + cos(2*M_PI*x(1))) / 4.0;
         }
         else
         {
            val = 100.0 + ((cos(2*M_PI*x(2)) + 2) *
                           (cos(2*M_PI*x(0)) + cos(2*M_PI*x(1))) - 2) / 16.0;
         }
         return val/denom;
      }
      case 1: return 0.0; // This case in initialized in main().
      case 2: if (x(0) < 0.5) { return 1.0 / rho0(x) / (gamma(x) - 1.0); }
         else { return 0.1 / rho0(x) / (gamma(x) - 1.0); }
      case 3: if (x(0) > 1.0) { return 0.1 / rho0(x) / (gamma(x) - 1.0); }
         else { return 1.0 / rho0(x) / (gamma(x) - 1.0); }
      case 4: return 1.0;
      case 5: return T_min + (T_max - T_min) *
                     (0.5 * (tanh(T_gradscale * (0.5 - x.Norml2())) + 1.0));
      case 6: return T_min + (T_max - T_min) *
                     (0.5 * (tanh(T_gradscale * (0.5 - x.Norml2())) + 1.0));
      default: MFEM_ABORT("Bad number given for problem id!"); return 0.0;
   }
}

} // namespace hydrodynamics

} // namespace mfem

#endif // MFEM_USE_MPI
