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

#include "m1_solver.hpp"

#ifdef MFEM_USE_MPI

using namespace std;

namespace mfem
{

namespace nth
{

void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    ParGridFunction &gf, const char *title,
                    int x, int y, int w, int h, bool vec)
{
   ParMesh &pmesh = *gf.ParFESpace()->GetParMesh();
   MPI_Comm comm = pmesh.GetComm();

   int num_procs, myid;
   MPI_Comm_size(comm, &num_procs);
   MPI_Comm_rank(comm, &myid);

   bool newly_opened = false;
   int connection_failed;

   do
   {
      if (myid == 0)
      {
         if (!sock.is_open() || !sock)
         {
            sock.open(vishost, visport);
            sock.precision(8);
            newly_opened = true;
         }
         sock << "solution\n";
      }

      pmesh.PrintAsOne(sock);
      gf.SaveAsOne(sock);

      if (myid == 0 && newly_opened)
      {
         sock << "window_title '" << title << "'\n"
              << "window_geometry "
              << x << " " << y << " " << w << " " << h << "\n"
              << "keys maaAcl";
         if ( vec ) { sock << "vvv"; }
         sock << endl;
      }

      if (myid == 0)
      {
         connection_failed = !sock && !newly_opened;
      }
      MPI_Bcast(&connection_failed, 1, MPI_INT, 0, comm);
   }
   while (connection_failed);
}

M1Operator::M1Operator(int size, 
                       ParFiniteElementSpace &h1_fes,
                       ParFiniteElementSpace &l2_fes, 
                       Array<int> &essential_tdofs,
                       ParGridFunction &rho0, 
                       double cfl_, 
                       NTHvHydroCoefficient *msp_,
                       NTHvHydroCoefficient *sourceI0_,
                       VectorCoefficient *Efield_,
					   VectorCoefficient *Bfield_,
					   ParGridFunction &x_gf_, 
                       ParGridFunction &T_gf_, 
                       bool pa, double cgt, int cgiter)
   : TimeDependentOperator(size),
     H1FESpace(h1_fes), L2FESpace(l2_fes),
     H1compFESpace(h1_fes.GetParMesh(), h1_fes.FEColl(), 1),
     ess_tdofs(essential_tdofs),
     dim(h1_fes.GetMesh()->Dimension()),
     nzones(h1_fes.GetMesh()->GetNE()),
     l2dofs_cnt(l2_fes.GetFE(0)->GetDof()),
     h1dofs_cnt(h1_fes.GetFE(0)->GetDof()),
     cfl(cfl_), p_assembly(pa), cg_rel_tol(cgt), cg_max_iter(cgiter),
     Mf1(&h1_fes), Mscattf1(&h1_fes), Bfieldf1(&h1_fes),
     MSf0(l2dofs_cnt, l2dofs_cnt, nzones),
     Mf0_inv(l2dofs_cnt, l2dofs_cnt, nzones),
     integ_rule(IntRules.Get(h1_fes.GetMesh()->GetElementBaseGeometry(),
                             3*h1_fes.GetOrder(0) + l2_fes.GetOrder(0) - 1)),
     quad_data(dim, nzones, integ_rule.GetNPoints()),
     quad_data_is_current(false),
     Divf0(&l2_fes, &h1_fes), Efieldf0(&l2_fes, &h1_fes),
     Divf1(&l2_fes, &h1_fes), AEfieldf1(&l2_fes, &h1_fes),
     AIEfieldf1(&l2_fes, &h1_fes),
     ForcePA(&quad_data, h1_fes, l2_fes),
     VMassPA(&quad_data, H1compFESpace), locEMassPA(&quad_data, l2_fes),
     locCG(), timer(), msp_pcf(msp_), sourceI0_pcf(sourceI0_), 
     Efield_pcf(Efield_), Bfield_pcf(Bfield_), x_gf(x_gf_)
{
   GridFunctionCoefficient rho_coeff(&rho0);

   // Values of rho0DetJ0 and Jac0inv at all quadrature points.
   const int nqp = integ_rule.GetNPoints();
   Vector rho_vals(nqp);
   for (int i = 0; i < nzones; i++)
   {
      rho0.GetValues(i, integ_rule, rho_vals);
      ElementTransformation *T = h1_fes.GetElementTransformation(i);
      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = integ_rule.IntPoint(q);
         T->SetIntPoint(&ip);

         DenseMatrixInverse Jinv(T->Jacobian());
         Jinv.GetInverseMatrix(quad_data.Jac0inv(i*nqp + q));

         const double rho0DetJ0 = T->Weight() * rho_vals(q);
         quad_data.rho0DetJ0w(i*nqp + q) = rho0DetJ0 *
                                           integ_rule.IntPoint(q).weight;
      }
   }

   // Initial local mesh size (assumes similar cells).
   double loc_area = 0.0, glob_area;
   int loc_z_cnt = nzones, glob_z_cnt;
   ParMesh *pm = H1FESpace.GetParMesh();
   for (int i = 0; i < nzones; i++) { loc_area += pm->GetElementVolume(i); }
   MPI_Allreduce(&loc_area, &glob_area, 1, MPI_DOUBLE, MPI_SUM, pm->GetComm());
   MPI_Allreduce(&loc_z_cnt, &glob_z_cnt, 1, MPI_INT, MPI_SUM, pm->GetComm());
   switch (pm->GetElementBaseGeometry(0))
   {
      case Geometry::SEGMENT:
         quad_data.h0 = glob_area / glob_z_cnt; break;
      case Geometry::SQUARE:
         quad_data.h0 = sqrt(glob_area / glob_z_cnt); break;
      case Geometry::TRIANGLE:
         quad_data.h0 = sqrt(2.0 * glob_area / glob_z_cnt); break;
      case Geometry::CUBE:
         quad_data.h0 = pow(glob_area / glob_z_cnt, 1.0/3.0); break;
      case Geometry::TETRAHEDRON:
         quad_data.h0 = pow(6.0 * glob_area / glob_z_cnt, 1.0/3.0); break;
      default: MFEM_ABORT("Unknown zone type!");
   }
   quad_data.h0 /= (double) H1FESpace.GetOrder(0);
 

   // Standard local assembly and inversion for energy mass matrices.
   DenseMatrix Mf0_(l2dofs_cnt);
   DenseMatrixInverse inv(&Mf0_);
   Mass0cIntegrator mi(quad_data);
   mi.SetIntRule(&integ_rule);
   for (int i = 0; i < nzones; i++)
   {
      mi.AssembleElementMatrix(*l2_fes.GetFE(i),
                               *l2_fes.GetElementTransformation(i), Mf0_);
      MSf0(i) = Mf0_;
      inv.Factor();
      inv.GetInverseMatrix(Mf0_inv(i));
   }

   // Standard assembly for the velocity mass matrix.
   Mass1NuIntegrator *f1mi = new Mass1NuIntegrator(quad_data);
   f1mi->SetIntRule(&integ_rule);
   Mf1.AddDomainIntegrator(f1mi);
   Mf1.Assemble();

   Mass1NutIntegrator *f1scati = new Mass1NutIntegrator(quad_data);
   f1scati->SetIntRule(&integ_rule);
   Mscattf1.AddDomainIntegrator(f1scati);
   Mscattf1.Assemble();

   BfieldIntegrator *f1bfi = new BfieldIntegrator(quad_data);
   f1bfi->SetIntRule(&integ_rule);
   Bfieldf1.AddDomainIntegrator(f1bfi);
   Bfieldf1.Assemble();

   Divf1Integrator *f1di = new Divf1Integrator(quad_data);
   f1di->SetIntRule(&integ_rule);
   Divf1.AddDomainIntegrator(f1di);
   // Make a dummy assembly to figure out the sparsity.
   Divf1.Assemble(0);
   Divf1.Finalize(0);

   Divf0Integrator *tfi = new Divf0Integrator(quad_data);
   tfi->SetIntRule(&integ_rule);
   Divf0.AddDomainIntegrator(tfi);
   // Make a dummy assembly to figure out the sparsity.
   Divf0.Assemble(0);
   Divf0.Finalize(0);

   EfieldIntegrator *f0ei = new EfieldIntegrator(quad_data);
   f0ei->SetIntRule(&integ_rule);
   Efieldf0.AddDomainIntegrator(f0ei);
   // Make a dummy assembly to figure out the sparsity.
   Efieldf0.Assemble(0);
   Efieldf0.Finalize(0);

   AEfieldIntegrator *f1aei = new AEfieldIntegrator(quad_data);
   f1aei->SetIntRule(&integ_rule);
   AEfieldf1.AddDomainIntegrator(f1aei);
   // Make a dummy assembly to figure out the sparsity.
   AEfieldf1.Assemble(0);
   AEfieldf1.Finalize(0);

   AIEfieldIntegrator *f1aiei = new AIEfieldIntegrator(quad_data);
   f1aiei->SetIntRule(&integ_rule);
   AIEfieldf1.AddDomainIntegrator(f1aiei);
   // Make a dummy assembly to figure out the sparsity.
   AIEfieldf1.Assemble(0);
   AIEfieldf1.Finalize(0);

   if (p_assembly)
   {
      if (tensors1D == NULL)
      {
	     tensors1D = new Tensors1D(H1FESpace.GetFE(0)->GetOrder(),
                                   L2FESpace.GetFE(0)->GetOrder(),
                                   int(floor(0.7 + pow(nqp, 1.0 / dim))));
	  }
      if (evaluator == NULL) { evaluator = new FastEvaluator(H1FESpace); }
   }

   locCG.SetOperator(locEMassPA);
   locCG.iterative_mode = false;
   locCG.SetRelTol(1e-8);
   locCG.SetAbsTol(1e-8 * numeric_limits<double>::epsilon());
   locCG.SetMaxIter(200);
   locCG.SetPrintLevel(0);
}

void M1Operator::Mult(const Vector &S, Vector &dS_dt) const
{
   dS_dt = 0.0;

   const double velocity = GetTime(); 

   UpdateQuadratureData(velocity, S);
   const double alphavT = msp_pcf->GetVelocityScale();
   const double velocity_scaled = velocity * alphavT;

   sourceI0_pcf->SetVelocity(velocity);
   ParGridFunction I0source(&L2FESpace);
   I0source.ProjectCoefficient(*sourceI0_pcf);

   // The monolithic BlockVector stores the unknown fields as follows:
   // - isotropic I0 (energy density)
   // - anisotropic I1 (flux density)

   const int VsizeL2 = L2FESpace.GetVSize();
   const int VsizeH1 = H1FESpace.GetVSize();

   ParGridFunction I0, I1;
   Vector* sptr = (Vector*) &S;
   I0.MakeRef(&L2FESpace, *sptr, 0);
   I1.MakeRef(&H1FESpace, *sptr, VsizeL2);

   ParGridFunction dI0, dI1;
   dI0.MakeRef(&L2FESpace, dS_dt, 0);
   dI1.MakeRef(&H1FESpace, dS_dt, VsizeL2);

   if (!p_assembly)
   {
      // Standard local assembly and inversion for energy mass matrices.
      DenseMatrix Mf0_(l2dofs_cnt);
      DenseMatrix explMf0_(l2dofs_cnt);
      DenseMatrixInverse inv(&explMf0_);
      ExplMass0Integrator explmi(quad_data);
      explmi.SetIntRule(&integ_rule);
      Mass0NuIntegrator mnui(quad_data);
      mnui.SetIntRule(&integ_rule);
      for (int i = 0; i < nzones; i++)
      {
         explmi.AssembleElementMatrix(*L2FESpace.GetFE(i),
                                      *L2FESpace.GetElementTransformation(i),
                                      explMf0_);
         inv.Factor();
         inv.GetInverseMatrix(Mf0_inv(i));
         mnui.AssembleElementMatrix(*L2FESpace.GetFE(i),
                                    *L2FESpace.GetElementTransformation(i),
                                    Mf0_);
         MSf0(i) = Mf0_;
      }

      Divf1 = 0.0;
      Divf0 = 0.0;
      AEfieldf1 = 0.0;
	  AIEfieldf1 = 0.0;
	  Efieldf0 = 0.0;
	  Mf1.Update();
      Bfieldf1.Update();
	  Mscattf1.Update();
      timer.sw_force.Start();
      Mf1.Assemble(); 
      Divf1.Assemble();
      Divf0.Assemble();
      AEfieldf1.Assemble(0);
	  AIEfieldf1.Assemble(0);
	  Efieldf0.Assemble(0);
	  Bfieldf1.Assemble();
	  Mscattf1.Assemble();
      timer.sw_force.Stop();
   }

   // Solve for df0dv.
   Array<int> l2dofs;
   Vector I0_rhs(VsizeL2), loc_rhs(l2dofs_cnt), loc_I0source(l2dofs_cnt),
          loc_MSf0MultI0source(l2dofs_cnt), loc_dI0(l2dofs_cnt);
   if (p_assembly)
   {
      timer.sw_force.Start();
      ForcePA.MultTranspose(I1, I0_rhs);
      timer.sw_force.Stop();
      timer.dof_tstep += L2FESpace.GlobalTrueVSize();

      for (int z = 0; z < nzones; z++)
      {
         L2FESpace.GetElementDofs(z, l2dofs);
         I0_rhs.GetSubVector(l2dofs, loc_rhs);
         locEMassPA.SetZoneId(z);
         //
         I0source.GetSubVector(l2dofs, loc_I0source);
         locEMassPA.Mult(loc_I0source, loc_MSf0MultI0source);
         loc_rhs += loc_MSf0MultI0source;
         //
         timer.sw_cgL2.Start();
         locCG.Mult(loc_rhs, loc_dI0);
         timer.sw_cgL2.Stop();
         timer.L2dof_iter += locCG.GetNumIterations() * l2dofs_cnt;
         dI0.SetSubVector(l2dofs, loc_dI0);
      }
   }
   else
   {
      timer.sw_force.Start();
      Divf0.MultTranspose(I1, I0_rhs);
      //Efieldf0.AddMultTranspose(I1, I0_rhs, 
      //                          2.0 / velocity_scaled / velocity_scaled);
      timer.sw_force.Stop();
      timer.dof_tstep += L2FESpace.GlobalTrueVSize();
      for (int z = 0; z < nzones; z++)
      {
         L2FESpace.GetElementDofs(z, l2dofs);
         I0_rhs.GetSubVector(l2dofs, loc_rhs);
         //
         I0source.GetSubVector(l2dofs, loc_I0source);
         MSf0(z).Mult(loc_I0source, loc_MSf0MultI0source);
         loc_rhs += loc_MSf0MultI0source;
         //
         timer.sw_cgL2.Start();
         // Scale rhs because of the normalized velocity, i.e. 
	     // Mf0*df0dv = 1/alphavT*Mf0*df0dvnorm = loc_rhs.
	     loc_rhs *= alphavT;
		 Mf0_inv(z).Mult(loc_rhs, loc_dI0);
         timer.sw_cgL2.Stop();
         timer.L2dof_iter += l2dofs_cnt;
         dI0.SetSubVector(l2dofs, loc_dI0);
		 //loc_dI0.Print();
      }
   }

   // Solve for df1dv.
   Vector rhs(VsizeH1), B, X;
   if (p_assembly)
   {
      timer.sw_force.Start();
      ForcePA.Mult(I0, rhs);
      timer.sw_force.Stop();
      timer.dof_tstep += H1FESpace.GlobalTrueVSize();
      rhs.Neg();

      // Partial assembly solve for each velocity component.
      const int size = H1compFESpace.GetVSize();
      for (int c = 0; c < dim; c++)
      {
         Vector rhs_c(rhs.GetData() + c*size, size),
                dI1_c(dI1.GetData() + c*size, size);

         Array<int> c_tdofs;
         Array<int> ess_bdr(H1FESpace.GetParMesh()->bdr_attributes.Max());
         // Attributes 1/2/3 correspond to fixed-x/y/z boundaries, i.e.,
         // we must enforce v_x/y/z = 0 for the velocity components.
         ess_bdr = 0; ess_bdr[c] = 1;
         // Essential true dofs as if there's only one component.
         H1compFESpace.GetEssentialTrueDofs(ess_bdr, c_tdofs);

         dI1_c = 0.0;
         Vector B(H1compFESpace.TrueVSize()), X(H1compFESpace.TrueVSize());
         H1compFESpace.Dof_TrueDof_Matrix()->MultTranspose(rhs_c, B);
         H1compFESpace.GetRestrictionMatrix()->Mult(dI1_c, X);

         VMassPA.EliminateRHS(c_tdofs, B);

         CGSolver cg(H1FESpace.GetParMesh()->GetComm());
         cg.SetOperator(VMassPA);
         cg.SetRelTol(cg_rel_tol);
         cg.SetAbsTol(0.0);
         cg.SetMaxIter(cg_max_iter);
         cg.SetPrintLevel(-1);
         timer.sw_cgH1.Start();
         cg.Mult(B, X);
         timer.sw_cgH1.Stop();
         timer.H1dof_iter += cg.GetNumIterations() *
                             H1compFESpace.GlobalTrueVSize();
         H1compFESpace.Dof_TrueDof_Matrix()->Mult(X, dI1_c);
      }
   }
   else
   {
      timer.sw_force.Start();
      Divf1.Mult(I0, rhs);
	  rhs.Neg();
      // dI0 negative (dfMdv) in diffusive regime. 
	  //AEfieldf1.AddMult(dI0, rhs, 1.0 / velocity_scaled);
      // Watch out! dI0 has been multiplied by alphavT because of it is 
	  // integrated along the normalized velocity dimension.
	  AEfieldf1.AddMult(dI0, rhs, 1.0 / velocity_scaled / alphavT);
	  //AEfieldf1.AddMult(I0source, rhs, 1.0 / velocity_scaled);
	  //AIEfieldf1.AddMult(I0, rhs, 1.0 / velocity_scaled / velocity_scaled);
	  Bfieldf1.AddMult(I1, rhs, 1.0 / velocity_scaled);
	  Mscattf1.AddMult(I1, rhs, 1.0 / velocity_scaled);
	  timer.sw_force.Stop();
      timer.dof_tstep += H1FESpace.GlobalTrueVSize();
      // Scale rhs because of the normalized velocity, i.e. 
	  // Mf1*df1dv = 1/alphavT*Mf1*df1dvnorm = rhs.
	  rhs *= alphavT;
      HypreParMatrix A;
      dI1 = 0.0;
      Mf1.FormLinearSystem(ess_tdofs, dI1, rhs, A, X, B);
      CGSolver cg(H1FESpace.GetParMesh()->GetComm());
      cg.SetOperator(A);
      cg.SetRelTol(1e-8); cg.SetAbsTol(0.0);
      cg.SetMaxIter(200);
      cg.SetPrintLevel(0);
      timer.sw_cgH1.Start();
      cg.Mult(B, X);
      timer.sw_cgH1.Stop();
      timer.H1dof_iter += cg.GetNumIterations() *
                          H1compFESpace.GlobalTrueVSize();
      Mf1.RecoverFEMSolution(X, rhs, dI1);
   }

   quad_data_is_current = false;
}

double M1Operator::GetVelocityStepEstimate(const Vector &S) const
{
   const double velocity = GetTime();
   UpdateQuadratureData(velocity, S);

   double glob_dt_est;
   MPI_Allreduce(&quad_data.dt_est, &glob_dt_est, 1, MPI_DOUBLE, MPI_MIN,
                 H1FESpace.GetParMesh()->GetComm());
   return glob_dt_est;
}

void M1Operator::ResetVelocityStepEstimate() const
{
   quad_data.dt_est = numeric_limits<double>::infinity();
}

void M1Operator::ComputeDensity(ParGridFunction &rho)
{
   rho.SetSpace(&L2FESpace);

   DenseMatrix Mrho(l2dofs_cnt);
   Vector rhs(l2dofs_cnt), rho_z(l2dofs_cnt);
   Array<int> dofs(l2dofs_cnt);
   DenseMatrixInverse inv(&Mrho);
   MassIntegrator mi(&integ_rule);
   DensityIntegrator di(quad_data);
   di.SetIntRule(&integ_rule);
   for (int i = 0; i < nzones; i++)
   {
      di.AssembleRHSElementVect(*L2FESpace.GetFE(i),
                                *L2FESpace.GetElementTransformation(i), rhs);
      mi.AssembleElementMatrix(*L2FESpace.GetFE(i),
                               *L2FESpace.GetElementTransformation(i), Mrho);
      inv.Factor();
      inv.Mult(rhs, rho_z);
      L2FESpace.GetElementDofs(i, dofs);
      rho.SetSubVector(dofs, rho_z);
   }
}

void M1Operator::PrintTimingData(bool IamRoot, int steps)
{
   double my_rt[5], rt_max[5];
   my_rt[0] = timer.sw_cgH1.RealTime();
   my_rt[1] = timer.sw_cgL2.RealTime();
   my_rt[2] = timer.sw_force.RealTime();
   my_rt[3] = timer.sw_qdata.RealTime();
   my_rt[4] = my_rt[0] + my_rt[2] + my_rt[3];
   MPI_Reduce(my_rt, rt_max, 5, MPI_DOUBLE, MPI_MAX, 0, H1FESpace.GetComm());

   double mydata[2], alldata[2];
   mydata[0] = timer.L2dof_iter;
   mydata[1] = timer.quad_tstep;
   MPI_Reduce(mydata, alldata, 2, MPI_DOUBLE, MPI_SUM, 0, H1FESpace.GetComm());

   if (IamRoot)
   {
      using namespace std;
      cout << endl;
      cout << "CG (H1) total time: " << rt_max[0] << endl;
      cout << "CG (H1) rate (megadofs x cg_iterations / second): "
           << 1e-6 * timer.H1dof_iter / rt_max[0] << endl;
      cout << endl;
      cout << "CG (L2) total time: " << rt_max[1] << endl;
      cout << "CG (L2) rate (megadofs x cg_iterations / second): "
           << 1e-6 * alldata[0] / rt_max[1] << endl;
      cout << endl;
      cout << "Divergences total time: " << rt_max[2] << endl;
      cout << "Divergences rate (megadofs x timesteps / second): "
           << 1e-6 * timer.dof_tstep / rt_max[2] << endl;
      cout << endl;
      cout << "UpdateQuadData total time: " << rt_max[3] << endl;
      cout << "UpdateQuadData rate (megaquads x timesteps / second): "
           << 1e-6 * alldata[1] / rt_max[3] << endl;
      cout << endl;
      cout << "Major kernels total time (seconds): " << rt_max[4] << endl;
      cout << "Major kernels total rate (megadofs x time steps / second): "
           << 1e-6 * H1FESpace.GlobalTrueVSize() * steps / rt_max[4] << endl;
   }
}

void M1Operator::UpdateQuadratureData(double velocity, const Vector &S) const
{
   if (quad_data_is_current) { return; }
   timer.sw_qdata.Start();

   msp_pcf->SetVelocity(velocity);
   const double alphavT = msp_pcf->GetVelocityScale();
   const double velocity_scaled = velocity * alphavT;
   const int nqp = integ_rule.GetNPoints();

   ParGridFunction I0, I1;
   Vector* sptr = (Vector*) &S;
   I0.MakeRef(&L2FESpace, *sptr, 0);
   I1.MakeRef(&H1FESpace, *sptr, L2FESpace.GetVSize());

   Vector vector_vals(h1dofs_cnt * dim);
   DenseMatrix Jpi(dim), Jinv(dim), I0stress(dim), I0stressJiT(dim),
               I1stress(dim), I1stressJiT(dim),
               vecvalMat(vector_vals.GetData(), h1dofs_cnt, dim);
   Array<int> L2dofs, H1dofs;

   // Batched computations are needed, because hydrodynamic codes usually
   // involve expensive computations of material properties. Although this
   // miniapp uses simple EOS equations, we still want to represent the batched
   // cycle structure.
   int nzones_batch = 3;
   const int nbatches =  nzones / nzones_batch + 1; // +1 for the remainder.
   int nqp_batch = nqp * nzones_batch;
   double *msp_b = new double[nqp_batch],
          *rho_b = new double[nqp_batch];
   // Electric and Magnetic fields for all quadrature points in the batch.
   //DenseMatrix *Efield_b = new DenseMatrix[nqp_batch],
   //            *Bfield_b = new DenseMatrix[nqp_batch];
   // Jacobians of reference->physical transformations for all quadrature
   // points in the batch.
   DenseTensor *Jpr_b = new DenseTensor[nqp_batch],
               *AM1_b = new DenseTensor[nqp_batch];
   for (int b = 0; b < nbatches; b++)
   {
      int z_id = b * nzones_batch; // Global index over zones.
      // The last batch might not be full.
      if (z_id == nzones) { break; }
      else if (z_id + nzones_batch > nzones)
      {
         nzones_batch = nzones - z_id;
         nqp_batch    = nqp * nzones_batch;
      }

      for (int z = 0; z < nzones_batch; z++)
      {
         ElementTransformation *T = H1FESpace.GetElementTransformation(z_id);
         //Efield_b[z].SetSize(dim, nqp);
		 //Bfield_b[z].SetSize(dim, nqp);
		 Jpr_b[z].SetSize(dim, dim, nqp);
         AM1_b[z].SetSize(dim, dim, nqp);

         if (p_assembly)
         {
            // All reference->physical Jacobians at the quadrature points.
            H1FESpace.GetElementVDofs(z_id, H1dofs);
            x_gf.GetSubVector(H1dofs, vector_vals);
            evaluator->GetVectorGrad(vecvalMat, Jpr_b[z]);
         }
         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = integ_rule.IntPoint(q);
            T->SetIntPoint(&ip);
            //Vector Efield;
            //Efield_pcf->Eval(Efield, *T, ip);
            //Efield_b[z].SetCol(q, Efield); 
            //Vector Bfield;
            //Bfield_pcf->Eval(Bfield, *T, ip);
            //Bfield_b[z].SetCol(q, Bfield);
            if (!p_assembly) { Jpr_b[z](q) = T->Jacobian(); }
            const double detJ = Jpr_b[z](q).Det();

            const int idx = z * nqp + q;
            //rho_b[idx] = msp_pcf->GetRho(*T, ip);
            rho_b[idx] = quad_data.rho0DetJ0w(z_id*nqp + q) / 
                         detJ / ip.weight;
            msp_b[idx] = msp_pcf->Eval(*T, ip, rho_b[idx]);
            // M1 closure.
            // Matric closure maximizing angular entropy
            // A = 1/3*I + M^2/2*(1 + M^2)*((f1xf1^T)/f1^2 - 1/3*I),
            // where M = |f1|/|f0| must be in (0, 1), "isotropic-freestreaming".
            // Isotropic matrix.
            DenseMatrix Iso;
            Iso.Diag(1.0 / 3.0, dim);
            AM1_b[z](q) = Iso;
/*
            double normlim = 1e-32;
            double anisolim = 1e-1;
            //double f0norm = I0.GetValue((*T).ElementNo, ip);
            double f0norm = abs(I0.GetValue((*T).ElementNo, ip));
            //double f0norm = max(normlim, I0.GetValue((*T).ElementNo, ip));
            Vector f1;
            I1.GetVectorValue((*T).ElementNo, ip, f1);
            double f1norm = f1.Norml2();
            if (f0norm < normlim || f1norm < normlim)
            {
               AM1_b[z](q) = Iso;
            }
            else if (f1norm / f0norm < anisolim)
            {
               AM1_b[z](q) = Iso;
            }
            else
            {
               double M = min(f1norm / f0norm, 1.0);
               double Msquare = M * M;
               double c = Msquare / 2.0 * (1.0 + Msquare);
               DenseMatrix normf1xf1T(dim);
               normf1xf1T.Diag(1.0, dim);
               //f1 *= 1.0 / f1norm;
			   // f1 directional matrix.
               //MultVVt(f1, normf1xf1T);
               //normf1xf1T = 0.0;
			   //for (int vd = 0; vd < dim; vd++)
               //{
               //   normf1xf1T(vd, vd) = f1(vd) * f1(vd);
               //}
               // Construct the closure matrix.
			   AM1_b[z](q) = 0.0;
               AM1_b[z](q).Add(1.0 - c, Iso);
               AM1_b[z](q).Add(c, normf1xf1T);
               //cout << "f0norm, f1norm, c:" << f0norm << ", " << f1norm
			   //     << ", " << c << endl << flush;
               //cout << "normf1xf1T:" << endl << flush;
               //normf1xf1T.Print();
               //cout << "AM1:" << endl << flush;
               //AM1_b[z](q).Print();
            }
*/
         }
         ++z_id;
      }

      z_id -= nzones_batch;
      for (int z = 0; z < nzones_batch; z++)
      {
         ElementTransformation *T = H1FESpace.GetElementTransformation(z_id);
         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = integ_rule.IntPoint(q);
            T->SetIntPoint(&ip);
            double f0min = 1e-32;
            double f0 = max(f0min, I0.GetValue((*T).ElementNo, ip));
            Vector f1;
            I1.GetVectorValue((*T).ElementNo, ip, f1);

            // Note that the Jacobian was already computed above. We've chosen
            // not to store the Jacobians for all batched quadrature points.
            const DenseMatrix &Jpr = Jpr_b[z](q);
            CalcInverse(Jpr, Jinv);
            const double detJ = Jpr.Det();
            const DenseMatrix &AM1 = AM1_b[z](q);
            DenseMatrix I; I.Diag(1.0, dim);
            const double msp = msp_b[z*nqp + q];
            double rho = rho_b[z*nqp + q];

            Vector Efield(dim), Bfield(dim), AEfield(dim), AIEfield(dim);
            // TODO Here the vector E and B evaluation.
            // And consequent evaluation of AE and AIE.
            //Efield_b[z].GetColumn(q, Efield);
            Efield_pcf->Eval(Efield, *T, ip);
			Bfield_pcf->Eval(Bfield, *T, ip);
			//Efield = 0.0;
            //Bfield = 0.0; 
            AM1.Mult(Efield, AEfield);
            AIEfield = 0.0;
            AM1.AddMult_a(3.0, Efield, AIEfield);
            I.AddMult_a(-1.0, Efield, AIEfield);
            // Time step estimate at the point. Here the more relevant length
            // scale is related to the actual mesh deformation; we use the min
            // singular value of the ref->physical Jacobian. In addition, the
            // time step estimate should be aware of the presence of shocks.
            const double h_min =
               Jpr.CalcSingularvalue(dim-1) / (double) H1FESpace.GetOrder(0);

            // The scaled cfl condition on velocity step.
            double dv = h_min * msp / alphavT; // / rho;
            quad_data.dt_est = min(quad_data.dt_est, cfl * dv);
            I0stress = 0.0;
            I1stress = 0.0;
			for (int d = 0; d < dim; d++)
            {
               I0stress(d, d) = 1.0;
               I1stress(d, d) = 1.0 / 3.0; // P1 closure.
            }
            // M1 closure. See the construction above.
            //I1stress = AM1;

            // Quadrature data for partial assembly of the force operator.
            MultABt(I0stress, Jinv, I0stressJiT);
            I0stressJiT *= integ_rule.IntPoint(q).weight * detJ;
            MultABt(I1stress, Jinv, I1stressJiT);
            I1stressJiT *= integ_rule.IntPoint(q).weight * detJ;
            for (int vd = 0 ; vd < dim; vd++)
            {
               for (int gd = 0; gd < dim; gd++)
               {
                  quad_data.stress1JinvT(vd)(z_id*nqp + q, gd) =
                     I1stressJiT(vd, gd);
                  quad_data.stress0JinvT(vd)(z_id*nqp + q, gd) =
                     I0stressJiT(vd, gd);
               }
               // Extensive vector quadrature data.
               quad_data.Einvrho(z_id*nqp + q, vd) = Efield(vd) / rho;
               quad_data.AEinvrho(z_id*nqp + q, vd) = AEfield(vd) / rho;
               //cout << "AE/rho: " << //Efield(0) 
			   //  quad_data.AEinvrho(z_id*nqp + q, vd) 
               //  << endl << flush;
               quad_data.AIEinvrho(z_id*nqp + q, vd) = AIEfield(vd) / rho;
               quad_data.Binvrho(z_id*nqp + q, vd) = Bfield(vd) / rho;
            }

            // Extensive scalar quadrature data.
            const double Zbar = 1.0;
            //const double Zbar = 10.0;
            //const double Zbar = 50.0;
            //const double Zbar = 100.0;
            // Introducing a AWBS correction factor in order to obtain 
            // an appropriate SH diffusive limit!
            const double corrAWBS = (688.9*Zbar + 114.4) / 
                                    (Zbar*Zbar + 1038.0*Zbar + 474.1);
            quad_data.nuinvrho(z_id*nqp + q) = corrAWBS * msp / rho; //nue/rho;
            //quad_data.nuinvrho(z_id*nqp + q) = msp / rho; //nue/rho;
            quad_data.Ef1invvf0rho(z_id*nqp + q) = Efield * f1
                                                   / velocity_scaled / f0
                                                   / rho;
            //cout << "Ef1/v/f0/rho: " <<  quad_data.Ef1invvf0rho(z_id*nqp + q) 
            //     << endl << flush;
            quad_data.nutinvrho(z_id*nqp + q) = Zbar * msp / rho;
         }
         ++z_id;
      }
   }
   delete [] msp_b;
   delete [] rho_b;
   delete [] Jpr_b;
   delete [] AM1_b;
   quad_data_is_current = true;

   timer.sw_qdata.Stop();
   timer.quad_tstep += nzones * nqp;
}

} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI
