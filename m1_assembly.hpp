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

#ifndef MFEM_M1_ASSEMBLY
#define MFEM_M1_ASSEMBLY

#include "mfem.hpp"


#ifdef MFEM_USE_MPI

#include <memory>
#include <iostream>

namespace mfem
{

namespace nth
{

// Container for all data needed at quadrature points.
struct QuadratureData
{
   // TODO: use QuadratureFunctions?

   // Reference to physical Jacobian for the initial mesh. These are computed
   // only at time zero and stored here.
   DenseTensor Jac0inv;

   // Quadrature data used for full/partial assembly of the force operator. At
   // each quadrature point, it combines the stress, inverse Jacobian,
   // determinant of the Jacobian and the integration weight. It must be
   // recomputed in every time step.
   DenseTensor stress1JinvT, stress0JinvT;

   // Quadrature data used for full/partial assembly of the mass matrices. At
   // time zero, we compute and store (rho0 * det(J0) * qp_weight) at each
   // quadrature point. Note the at any other time, we can compute
   // rho = rho0 * det(J0) / det(J), representing the notion of pointwise mass
   // conservation.
   Vector rho0DetJ0w;

   // The pointwise equality rho * detJ = rho0 * detJ0 is used by integrators.
   // Electric and magnetic fields. 
   DenseMatrix Einvvnue, AEinvvnue, AIEinvv2nue, Binvvnue;
   // Angular scattering integrator.
   Vector nutinvvnue;
   // Explicit zero moment "mass" integrator.
   Vector Ef1invvnuef0;
   
   // Complementary matrix to divergence. 
   DenseMatrix invrhoAgradrhoinvnue;

   // Initial length scale. This represents a notion of local mesh size. We
   // assume that all initial zones have similar size.
   double h0;

   // Estimate of the minimum time step over all quadrature points. This is
   // recomputed at every time step to achieve adaptive time stepping.
   double dt_est;

   QuadratureData(int dim, int nzones, int quads_per_zone)
      : Jac0inv(dim, dim, nzones * quads_per_zone),
        stress1JinvT(nzones * quads_per_zone, dim, dim),
        stress0JinvT(nzones * quads_per_zone, dim, dim),
		rho0DetJ0w(nzones * quads_per_zone),
        Einvvnue(nzones * quads_per_zone, dim),
        AEinvvnue(nzones * quads_per_zone, dim),
        AIEinvv2nue(nzones * quads_per_zone, dim),
        Binvvnue(nzones * quads_per_zone, dim),
        nutinvvnue(nzones * quads_per_zone),
        Ef1invvnuef0(nzones * quads_per_zone),
        invrhoAgradrhoinvnue(nzones * quads_per_zone, dim) { }
};

// Stores values of the one-dimensional shape functions and gradients at all 1D
// quadrature points. All sizes are (dofs1D_cnt x quads1D_cnt).
struct Tensors1D
{
   // H1 shape functions and gradients, L2 shape functions.
   DenseMatrix HQshape1D, HQgrad1D, LQshape1D;

   Tensors1D(int H1order, int L2order, int nqp1D);
};
extern const Tensors1D *tensors1D;

class FastEvaluator
{
   const int dim;
   ParFiniteElementSpace &H1FESpace;

public:
   FastEvaluator(ParFiniteElementSpace &h1fes)
      : dim(h1fes.GetMesh()->Dimension()), H1FESpace(h1fes) { }

   void GetL2Values(const Vector &vecL2, Vector &vecQP) const;
   // The input vec is an H1 function with dim components, over a zone.
   // The output is J_ij = d(vec_i) / d(x_j) with ij = 1 .. dim.
   void GetVectorGrad(const DenseMatrix &vec, DenseTensor &J) const;
};
extern const FastEvaluator *evaluator;

// This class is used only for visualization. It assembles (rho, phi) in each
// zone, which is used by LagrangianHydroOperator::ComputeDensity to do an L2
// projection of the density.
class DensityIntegrator : public LinearFormIntegrator
{
private:
   const QuadratureData &quad_data;

public:
   DensityIntegrator(QuadratureData &quad_data_) : quad_data(quad_data_) { }

   virtual void AssembleRHSElementVect(const FiniteElement &fe,
                                       ElementTransformation &Tr,
                                       Vector &elvect);
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class M0Integrator : public BilinearFormIntegrator
{
protected:
   const QuadratureData &quad_data;

public:
   M0Integrator(QuadratureData &quad_data_) : quad_data(quad_data_) { }

   virtual void AssembleElementMatrix(const FiniteElement &fe,
                                      ElementTransformation &Trans,
                                      DenseMatrix &elmat)
   { AssembleElementMatrix2(fe, fe, Trans, elmat); }
   virtual void AssembleElementMatrix2(const FiniteElement &trial_fe,
                                       const FiniteElement &test_fe,
                                       ElementTransformation &Trans,
                                       DenseMatrix &elmat);
   virtual double GetIntegrator(int i) = 0;
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class Mass1Integrator : public BilinearFormIntegrator
{
protected:
   const QuadratureData &quad_data;

public:
   Mass1Integrator(QuadratureData &quad_data_) : quad_data(quad_data_) { }

   virtual void AssembleElementMatrix(const FiniteElement &fe,
                                      ElementTransformation &Trans,
                                      DenseMatrix &elmat)
   { AssembleElementMatrix2(fe, fe, Trans, elmat); }
   virtual void AssembleElementMatrix2(const FiniteElement &trial_fe,
                                       const FiniteElement &test_fe,
                                       ElementTransformation &Trans,
                                       DenseMatrix &elmat);
   virtual double GetIntegrator(int i) = 0;
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class DivIntegrator : public BilinearFormIntegrator
{
protected:
   const QuadratureData &quad_data;

public:
   DivIntegrator(QuadratureData &quad_data_) : quad_data(quad_data_) { }

   virtual void AssembleElementMatrix2(const FiniteElement &trial_fe,
                                       const FiniteElement &test_fe,
                                       ElementTransformation &Trans,
                                       DenseMatrix &elmat);
   virtual double GetIntegrator(int q, int vd, int gd) = 0;
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class VdotIntegrator : public BilinearFormIntegrator
{
protected:
   const QuadratureData &quad_data;

public:
   VdotIntegrator(QuadratureData &quad_data_) : quad_data(quad_data_) { }

   virtual void AssembleElementMatrix2(const FiniteElement &trial_fe,
                                       const FiniteElement &test_fe,
                                       ElementTransformation &Trans,
                                       DenseMatrix &elmat);
   virtual double GetIntegrator(int q, int vd) = 0;
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class VcrossIntegrator : public BilinearFormIntegrator
{
protected:
   const QuadratureData &quad_data;

public:
   VcrossIntegrator(QuadratureData &quad_data_) : quad_data(quad_data_) { }

   virtual void AssembleElementMatrix(const FiniteElement &fe,
                                      ElementTransformation &Trans,
                                      DenseMatrix &elmat)
   { AssembleElementMatrix2(fe, fe, Trans, elmat); }
   virtual void AssembleElementMatrix2(const FiniteElement &trial_fe,
                                       const FiniteElement &test_fe,
                                       ElementTransformation &Trans,
                                       DenseMatrix &elmat);
   virtual double GetIntegrator(int q, int vd) = 0;
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class M0cIntegrator : public M0Integrator
{
private:
public:
   M0cIntegrator(QuadratureData &quad_data_) : M0Integrator(quad_data_) { }

   double GetIntegrator(int i);
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class ExplM0Integrator : public M0Integrator
{
private:
public:
   ExplM0Integrator(QuadratureData &quad_data_) : M0Integrator(quad_data_) { }

   double GetIntegrator(int i);
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class Mass1cIntegrator : public Mass1Integrator
{
private:
public:
   Mass1cIntegrator(QuadratureData &quad_data_) :
      Mass1Integrator(quad_data_) { }

   double GetIntegrator(int i);
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class Mass1NuIntegrator : public Mass1Integrator
{
private:
public:
   Mass1NuIntegrator(QuadratureData &quad_data_) :
      Mass1Integrator(quad_data_) { }

   double GetIntegrator(int i);
};

// Assembles element contributions to the global temperature force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class Divf0Integrator : public DivIntegrator
{
private:

public:
   Divf0Integrator(QuadratureData &quad_data_) : DivIntegrator(quad_data_) { }

   double GetIntegrator(int q, int vd, int gd);
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class Divf1Integrator : public DivIntegrator
{
private:
public:
   Divf1Integrator(QuadratureData &quad_data_) : DivIntegrator(quad_data_) { }

   double GetIntegrator(int q, int vd, int gd);
};

// Assembles element contributions to the global temperature force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class AgradIntegrator : public VdotIntegrator
{
private:

public:
   AgradIntegrator(QuadratureData &quad_data_) : VdotIntegrator(quad_data_) { }

   double GetIntegrator(int q, int vd);
};

// Assembles element contributions to the global temperature force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class EfieldIntegrator : public VdotIntegrator
{
private:

public:
   EfieldIntegrator(QuadratureData &quad_data_) : VdotIntegrator(quad_data_) { }

   double GetIntegrator(int q, int vd);
};

// Assembles element contributions to the global temperature force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class AEfieldIntegrator : public VdotIntegrator
{
private:

public:
   AEfieldIntegrator(QuadratureData &quad_data_) :
      VdotIntegrator(quad_data_) { }

   double GetIntegrator(int q, int vd);
};

// Assembles element contributions to the global temperature force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class AIEfieldIntegrator : public VdotIntegrator
{
private:

public:
   AIEfieldIntegrator(QuadratureData &quad_data_) :
      VdotIntegrator(quad_data_) { }

   double GetIntegrator(int q, int vd);
};

// Assembles element contributions to the global temperature force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class BfieldIntegrator : public VcrossIntegrator
{
private:

public:
   BfieldIntegrator(QuadratureData &quad_data_) :
      VcrossIntegrator(quad_data_) { }

   double GetIntegrator(int q, int vd);
};

// Performs partial assembly, which corresponds to (and replaces) the use of the
// LagrangianHydroOperator::Force global matrix.
class ForcePAOperator : public Operator
{
private:
   const int dim, nzones;

   QuadratureData *quad_data;
   ParFiniteElementSpace &H1FESpace, &L2FESpace;

   // Force matrix action on quadrilateral elements in 2D
   void MultQuad(const Vector &vecL2, Vector &vecH1) const;
   // Force matrix action on hexahedral elements in 3D
   void MultHex(const Vector &vecL2, Vector &vecH1) const;

   // Transpose force matrix action on quadrilateral elements in 2D
   void MultTransposeQuad(const Vector &vecH1, Vector &vecL2) const;
   // Transpose force matrix action on hexahedral elements in 3D
   void MultTransposeHex(const Vector &vecH1, Vector &vecL2) const;

public:
   ForcePAOperator(QuadratureData *quad_data_,
                   ParFiniteElementSpace &h1fes, ParFiniteElementSpace &l2fes)
      : dim(h1fes.GetMesh()->Dimension()), nzones(h1fes.GetMesh()->GetNE()),
        quad_data(quad_data_), H1FESpace(h1fes), L2FESpace(l2fes) { }

   virtual void Mult(const Vector &vecL2, Vector &vecH1) const;
   virtual void MultTranspose(const Vector &vecH1, Vector &vecL2) const;

   ~ForcePAOperator() { }
};

// Performs partial assembly for the velocity mass matrix.
class MassPAOperator : public Operator
{
private:
   const int dim, nzones;

   QuadratureData *quad_data;
   ParFiniteElementSpace &FESpace;

   Array<int> *ess_tdofs;

   mutable ParGridFunction x_gf, y_gf;

   // Mass matrix action on quadrilateral elements in 2D.
   void MultQuad(const Vector &x, Vector &y) const;
   // Mass matrix action on hexahedral elements in 3D.
   void MultHex(const Vector &x, Vector &y) const;

public:
   MassPAOperator(QuadratureData *quad_data_, ParFiniteElementSpace &fes)
      : Operator(fes.TrueVSize()),
        dim(fes.GetMesh()->Dimension()), nzones(fes.GetMesh()->GetNE()),
        quad_data(quad_data_), FESpace(fes), ess_tdofs(NULL),
        x_gf(&fes), y_gf(&fes)
   { }

   // Mass matrix action. We work with one velocity component at a time.
   virtual void Mult(const Vector &x, Vector &y) const;

   void EliminateRHS(Array<int> &dofs, Vector &b)
   {
      ess_tdofs = &dofs;
      for (int i = 0; i < dofs.Size(); i++) { b(dofs[i]) = 0.0; }
   }
};

// Performs partial assembly for the energy mass matrix on a single zone.
// Used to perform local CG solves, thus avoiding unnecessary communication.
class LocalMassPAOperator : public Operator
{
private:
   const int dim;
   int zone_id;

   QuadratureData *quad_data;

   // Mass matrix action on a quadrilateral element in 2D.
   void MultQuad(const Vector &x, Vector &y) const;
   // Mass matrix action on a hexahedral element in 3D.
   void MultHex(const Vector &x, Vector &y) const;

public:
   LocalMassPAOperator(QuadratureData *quad_data_, ParFiniteElementSpace &fes)
      : Operator(fes.GetFE(0)->GetDof()),
        dim(fes.GetMesh()->Dimension()), zone_id(0),
        quad_data(quad_data_)
   { }
   void SetZoneId(int zid) { zone_id = zid; }

   virtual void Mult(const Vector &x, Vector &y) const;
};

} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_M1_ASSEMBLY
