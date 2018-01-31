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

#ifndef MFEM_NTH_EOS
#define MFEM_NTH_EOS

#include "mfem.hpp"
#include "physics.hpp"

#ifdef MFEM_USE_MPI

namespace mfem
{

namespace nth
{

// Generic hydro equation of state (EOS) class,
// providing any physics related evaluation needed in NTH.
class IGEOS : public EOS
{
protected:
   // Mean ionization, even though a constant.
   double Zbar;
public:
   IGEOS(double kB_, double me_) : EOS(me_, kB_, 1.0, 1.0, 1.0) { Zbar = 1.0; }
   // Get Thermodynamic values.
   virtual double GetZbar(double index, double rho, double Te) { return Zbar; }
   virtual double GetPe(double index, double rho, double Te) {}
   virtual double GetPi(double index, double rho, double Te) {}
   virtual double GetEe(double index, double rho, double Te) {}
   virtual double GetEi(double index, double rho, double Te) {}
   virtual double GetdEedTe(double index, double rho, double Te) {}
   virtual double GetdEidTi(double index, double rho, double Te) {}
   virtual double GetdEedrho(double index, double rho, double Te) {}
   virtual double GetdEidrho(double index, double rho, double Te) {}
   virtual double GetSoundSpeed(double index, double rho, double Te) {}
   // IG specific functions.
   void SetZbar(double Zbar_) { Zbar = Zbar_; }
};

} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_NTH_EOS
