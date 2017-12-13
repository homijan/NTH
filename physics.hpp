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

#ifndef MFEM_NTH_PHYSICS
#define MFEM_NTH_PHYSICS

#include "mfem.hpp"

#ifdef MFEM_USE_MPI

//#include <memory>
//#include <iostream>
//#include <fstream>

namespace mfem
{

namespace nth
{

// Generic hydro equation of state (EOS) class,
// providing any physics related evaluation needed in NTH.
class EOS
{
protected:
   // Fundamental constants of nature.
   double kB, c, hbar, G;
   // Corresponding masses of electron and proton.
   double me;
public:
   EOS(double me_ = 1.0, double kB_ = 1.0, double c_ = 1.0,
       double hbar_ = 1.0, double G_ = 1.0)
      { kB = kB_, c = c_, hbar = hbar_, G = G_; me = me_; }
   // Get fundamental physical constants.
   double GetkB() {return kB; }
   double Getc() {return c; }
   double Getkhbar() {return hbar; }
   double GetG() {return G; }
   double Getme() {return me; }
   // Thermal velocity.
   double GetvTe(double Te) { return sqrt(kB * Te / me); }
   // Get Thermodynamic values.
   virtual double GetPe(double index, double rho, double Te) = 0;
   virtual double GetPi(double index, double rho, double Te) = 0;
   virtual double GetEe(double index, double rho, double Te) = 0;
   virtual double GetEi(double index, double rho, double Te) = 0;
   virtual double GetdEedTe(double index, double rho, double Te) = 0;
   virtual double GetdEidTi(double index, double rho, double Te) = 0;
   virtual double GetdEedrho(double index, double rho, double Te) = 0;
   virtual double GetdEidrho(double index, double rho, double Te) = 0;
   virtual double GetSoundSpeed(double index, double rho, double Te) = 0;

   ~EOS() { }
};

// Generic hydro coefficient.
class HydroCoefficient : public Coefficient
{
protected:
   // Fluid quantities used in calculations of physics.
   ParGridFunction &rho_gf, &Te_gf, &v_gf;
   // Space dependent material coefficient.
   Coefficient *material_pcf;
   // General equation of state.
   EOS *eos;
public:
   HydroCoefficient(ParGridFunction &rho_, ParGridFunction &Te_,
                    ParGridFunction &v_, Coefficient *material_, EOS *eos_)
      : rho_gf(rho_), Te_gf(Te_), v_gf(v_), material_pcf(material_), eos(eos_)
      { }
   virtual double Eval(ElementTransformation &T,
      const IntegrationPoint &ip) = 0;

   virtual ~HydroCoefficient() {};
};

// NTH hydro coefficient including velocity dependence.
class NTHvHydroCoefficient : public HydroCoefficient
{
   void SetVelocityScale(double alpha_, double Tmax)
      { alphavT = alpha * eos->GetvTe(Tmax); }
protected:
   // Velocity is always scaled wit respect to maximum thermal velocity
   // (its multiple) so it is in (0, 1)
   double alpha, Tmax, alphavT;
   // Current particle velocity from the velocity spectra.
   double velocity;
public:
   NTHvHydroCoefficient(ParGridFunction &rho_, ParGridFunction &T_,
                        ParGridFunction &v_, Coefficient *material_, EOS *eos_)
      : HydroCoefficient(rho_, T_, v_, material_, eos_)
	  { alpha = 1.0; Tmax = 1.0; SetVelocityScale(alpha, Tmax); }
   virtual double Eval(ElementTransformation &T,
      const IntegrationPoint &ip) = 0;
   virtual double Eval(ElementTransformation &T,
      const IntegrationPoint &ip, double rho) = 0;
   void SetVelocity(double v_) { velocity = v_; }
   void SetThermalVelocityMultiple(double alpha_)
      { alpha = alpha_; SetVelocityScale(alpha, Tmax); }
   void SetTmax(double Tmax_)
      { Tmax = Tmax_; SetVelocityScale(alpha, Tmax); }
   double GetVelocityScale() { return alphavT; }
   double GetRho(ElementTransformation &T, const IntegrationPoint &ip)
      { return rho_gf.GetValue(T.ElementNo, ip); }
};

// Classical mean-stopping-power coefficient.
class ClassicalMeanStoppingPower : public NTHvHydroCoefficient
{
protected:
public:
   ClassicalMeanStoppingPower(ParGridFunction &rho_, ParGridFunction &Te_,
                              ParGridFunction &v_, Coefficient *material_,
                              EOS *eos_)
      : NTHvHydroCoefficient(rho_, Te_, v_, material_, eos_) {}
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip,
                       double rho);
};

// General mean-free-path coefficient.
class MeanFreePath
{
protected:
public:
   MeanFreePath() {}
   virtual double EvalThermalMFP(ElementTransformation &T, 
                                 const IntegrationPoint &ip) = 0;
};

// Classical mean-free-path coefficient.
class ClassicalMeanFreePath : public MeanFreePath, 
                              public ClassicalMeanStoppingPower
{
protected:
public:
   ClassicalMeanFreePath(ParGridFunction &rho_, ParGridFunction &Te_,
                         ParGridFunction &v_, Coefficient *material_, 
                         EOS *eos_)
      : ClassicalMeanStoppingPower(rho_, Te_, v_, material_, eos_) {}
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip,
                       double rho);
   virtual double EvalThermalMFP(ElementTransformation &T, 
                                 const IntegrationPoint &ip);
};

// Classical Kn(mean-stopping-power) coefficient.
class KnudsenNumber : public HydroCoefficient
{
protected:
   MeanFreePath *mfp;
public:
   KnudsenNumber(ParGridFunction &rho_, ParGridFunction &Te_,
                 ParGridFunction &v_, Coefficient *material_,
                 EOS *eos_, MeanFreePath *mfp_)
      : HydroCoefficient(rho_, Te_, v_, material_, eos_), mfp(mfp_) {}
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
};

// AWBS source coefficient.
class AWBSI0Source : public NTHvHydroCoefficient
{
protected:
public:
   AWBSI0Source(ParGridFunction &rho_, ParGridFunction &Te_,
                ParGridFunction &v_, Coefficient *material_, EOS *eos_)
      : NTHvHydroCoefficient(rho_, Te_, v_, material_, eos_) {}
   double Eval(ElementTransformation &T, const IntegrationPoint &ip);
   double Eval(ElementTransformation &T, const IntegrationPoint &ip,
               double rho);
};

extern double a0;

} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_NTH_PHYSICS
