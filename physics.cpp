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
   double velocity_real = alphavT * velocity;
   double nu = a0 * rho / pow(velocity_real, 3.0);

   return nu;
}

double ClassicalMeanFreePath::EvalThermalMFP(ElementTransformation &T,
                                             const IntegrationPoint &ip,
                                             double rho)
{
   //double rho = rho_gf.GetValue(T.ElementNo, ip);
   double Te = Te_gf.GetValue(T.ElementNo, ip);
   // Set the scaled velocity to correspond to the local thermal velocity.
   velocity = eos->GetvTe(Te) / alphavT;
   double velocity_real = alphavT * velocity;
   // Compute the mean free path.
   double nu = ClassicalMeanStoppingPower::Eval(T, ip, rho);
   double mfp = velocity_real / nu;

   return mfp;
}

double KnudsenNumber::Eval(ElementTransformation &T,
                           const IntegrationPoint &ip)
{
   double rho = rho_gf.GetValue(T.ElementNo, ip);
   double Te = Te_gf.GetValue(T.ElementNo, ip);
   // Compute the mean free path.
   double lambda = mfp->EvalThermalMFP(T, ip, rho);
   // Compute the temperature and density length scales.
   Vector grad_Te;
   Te_gf.GetGradient(T, grad_Te);
   double L_Te = Te / grad_Te.Norml2();
   Vector grad_rho;
   rho_gf.GetGradient(T, grad_rho);
   double L_rho = rho / grad_rho.Norml2();
   // Return the Knudsen number of thermal velocity particle.
   return lambda / min(L_Te, L_rho);
}

void LorentzEfield::Eval(Vector &V, ElementTransformation &T,
                           const IntegrationPoint &ip)
{
   double rho = rho_gf.GetValue(T.ElementNo, ip);
   //double Te = Te_gf.GetValue(T.ElementNo, ip);
   double Te = max(1e-10, Te_gf.GetValue(T.ElementNo, ip));
   // Set the scaled velocity to correspond to the local thermal velocity.
   double vTe = eos->GetvTe(Te);
   // Compute the temperature and density length scales.
   Vector grad_Te;
   Te_gf.GetGradient(T, grad_Te);
   Vector grad_rho;
   rho_gf.GetGradient(T, grad_rho);
   // Return the Lorentz quasi-neutral (zero current) Efield.
   //Te = 1e4;
   //grad_Te = -1.0;
   V = grad_Te;
   V *= 2.5 / Te;
   //V *= rho;
   //V += grad_rho;
   //V *= 1.0 / rho;
   V *= pow(vTe, 2.0);
   //V *= 1e-16 * vTe * vTe;
   //V *= 4.55e-3 * v_th * v_th; // Reference -p 5 (T in (1000, 100))
   //V *= 6.085e-3 * v_th * v_th; // Reference -p 8 (T=5e2)
   //V *= 2.02e-3 * v_th * v_th; // Reference -p 8 (T=5e3)
   //V *= 1.4284e-3 * v_th * v_th; // Reference -p 8 (T=1e4)
   //V *= 1.01e-3 * v_th * v_th; // Reference -p 8 2xT (T=2e4)
   //V *= 0.71425e-3 * v_th * v_th; // Reference -p 8 4xT (T=4e4)
   //V *= 0.5051e-3 * v_th * v_th; // Reference -p 8 8xT (T=8e4)
   //V *= 0.35715e-3 * v_th * v_th; // Reference -p 8 16xT (T=16e4)
   //V *= 0.2525e-3 * v_th * v_th; // Reference -p 8 32xT (T=32e4)
   //V *= 0.127e-3 * v_th * v_th; // Reference -p 8 128xT (T=128e4)
}

double LorentzEfield::Eval(ElementTransformation &T,
                           const IntegrationPoint &ip)
{
   // Return the Lorentz quasi-neutral (zero current) Efield.
   Vector V;
   Eval(V, T, ip);

   return V.Norml2();
}

double AWBSI0Source::Eval(ElementTransformation &T,
                          const IntegrationPoint &ip, double rho)
{
   double pi = 3.14159265359;
   double Te = max(1e-10, Te_gf.GetValue(T.ElementNo, ip));
   double vTe = eos->GetvTe(Te);
   double velocity_real = alphavT * velocity;

   // Maxwell-Boltzmann distribution fM = ne*vT^3*(2/pi)^1.5*exp(-v^2/2/vT^2)
   double fM = 4.0 * pi *
               rho / pow(vTe, 3.0) / pow(2.0 * pi, 1.5) *
               exp(- pow(velocity_real, 2.0) / 2.0 / pow(vTe, 2.0));
   double dfMdv = - velocity_real / pow(vTe, 2.0) * fM;

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
