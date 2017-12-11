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

#include "physics.hpp"

#ifdef MFEM_USE_MPI

using namespace std;

namespace mfem
{

namespace nth
{

double a0 = 5e3;

double ClassicalMeanStoppingPower::Eval(ElementTransformation &T,
                                        const IntegrationPoint &ip, double rho)
{
   double Te = Te_gf.GetValue(T.ElementNo, ip);
   //double a = a0 * (Tmax * Tmax); //1e8; // The plasma collision model.
   double nu = a0 * rho / (pow(alphavT, 3.0) * pow(velocity, 3.0));

   return nu;
}

double ClassicalMeanStoppingPower::Eval(ElementTransformation &T,
                                        const IntegrationPoint &ip)
{
   double rho = rho_gf.GetValue(T.ElementNo, ip);

   return Eval(T, ip, rho);
}

double AWBSI0Source::Eval(ElementTransformation &T,
                          const IntegrationPoint &ip, double rho)
{
   double pi = 3.14159265359;
   double Te = max(1e-10, Te_gf.GetValue(T.ElementNo, ip));
   double vTe = eos->GetvTe(Te);

   // Maxwell-Boltzmann distribution fM = ne*vT^3*(2/pi)^1.5*exp(-v^2/2/vT^2)
   double fM = rho / pow(vTe, 3.0) / pow(2.0 * pi, 1.5) *
               exp(- pow(alphavT, 2.0) / 2.0 / pow(vTe, 2.0) *
               pow(velocity, 2.0));
   double dfMdv = - alphavT * velocity / pow(vTe, 2.0) * fM;

   return dfMdv;
}

double AWBSI0Source::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
   double rho = rho_gf.GetValue(T.ElementNo, ip);

   return Eval(T, ip, rho);
}


} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI
