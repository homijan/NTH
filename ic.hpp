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

#ifndef MFEM_NTH_IC
#define MFEM_NTH_IC

#include "mfem.hpp"

#ifdef MFEM_USE_MPI

namespace mfem
{

namespace nth
{

extern int nth_problem;

} // namespace nth

namespace hydrodynamics
{

extern int hydro_problem;
extern double T_max, T_min;
extern double rho_max, rho_min;
extern double T_gradscale, rho_gradscale;

double rho0(const Vector &x);
double gamma(const Vector &x);
void v0(const Vector &x, Vector &v);
double e0(const Vector &x);

} // namespace hydrodynamics

} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_NTH_EOS
