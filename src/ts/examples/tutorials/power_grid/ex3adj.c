
static char help[] = "Basic equation for generator stability analysis.\n";

/*F

\begin{eqnarray}
                 \frac{d \theta}{dt} = \omega_b (\omega - \omega_s)
                 \frac{2 H}{\omega_s}\frac{d \omega}{dt} & = & P_m - P_max \sin(\theta) -D(\omega - \omega_s)\\
\end{eqnarray}



  Ensemble of initial conditions
   ./ex2 -ensemble -ts_monitor_draw_solution_phase -1,-3,3,3      -ts_adapt_dt_max .01  -ts_monitor -ts_type rosw -pc_type lu -ksp_type preonly

  Fault at .1 seconds
   ./ex2           -ts_monitor_draw_solution_phase .42,.95,.6,1.05 -ts_adapt_dt_max .01  -ts_monitor -ts_type rosw -pc_type lu -ksp_type preonly

  Initial conditions same as when fault is ended
   ./ex2 -u 0.496792,1.00932 -ts_monitor_draw_solution_phase .42,.95,.6,1.05  -ts_adapt_dt_max .01  -ts_monitor -ts_type rosw -pc_type lu -ksp_type preonly 


F*/

/*
   Include "petscts.h" so that we can use TS solvers.  Note that this
   file automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
     petscksp.h   - linear solvers
*/
#include <petscts.h>

typedef struct {
  PetscScalar H,D,omega_b,omega_s,Pmax,Pm,E,V,X,u_s,c;
  PetscInt    beta;
  PetscReal   tf,tcl;
} AppCtx;

#undef __FUNCT__
#define __FUNCT__ "IFunction"
/*
     Defines the ODE passed to the ODE solver
*/
static PetscErrorCode IFunction(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,AppCtx *ctx)
{
  PetscErrorCode ierr;
  PetscScalar    *u,*udot,*f,Pmax;

  PetscFunctionBegin;
  /*  The next three lines allow us to access the entries of the vectors directly */
  ierr = VecGetArray(U,&u);CHKERRQ(ierr);
  ierr = VecGetArray(Udot,&udot);CHKERRQ(ierr);
  ierr = VecGetArray(F,&f);CHKERRQ(ierr);
  if ((t > ctx->tf) && (t < ctx->tcl)) Pmax = 0.0; /* A short-circuit on the generator terminal that drives the electrical power output (Pmax*sin(delta)) to 0 */
  else Pmax = ctx->Pmax;
  
  f[0] = udot[0] - ctx->omega_b*(u[1] - ctx->omega_s);
  f[1] = 2.0*ctx->H/ctx->omega_s*udot[1] +  Pmax*PetscSinScalar(u[0]) + ctx->D*(u[1] - ctx->omega_s)- ctx->Pm;

  ierr = VecRestoreArray(U,&u);CHKERRQ(ierr);
  ierr = VecRestoreArray(Udot,&udot);CHKERRQ(ierr);
  ierr = VecRestoreArray(F,&f);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IJacobian"
/*
     Defines the Jacobian of the ODE passed to the ODE solver. See TSSetIJacobian() for the meaning of a and the Jacobian.
*/
static PetscErrorCode IJacobian(TS ts,PetscReal t,Vec U,Vec Udot,PetscReal a,Mat A,Mat B,AppCtx *ctx)
{
  PetscErrorCode ierr;
  PetscInt       rowcol[] = {0,1};
  PetscScalar    *u,*udot,J[2][2],Pmax;

  PetscFunctionBegin;
  ierr = VecGetArray(U,&u);CHKERRQ(ierr);
  ierr = VecGetArray(Udot,&udot);CHKERRQ(ierr);
  if ((t > ctx->tf) && (t < ctx->tcl)) Pmax = 0.0; /* A short-circuit on the generator terminal that drives the electrical power output (Pmax*sin(delta)) to 0 */
  else Pmax = ctx->Pmax;

  J[0][0] = a;                       J[0][1] = -ctx->omega_b;
  J[1][1] = 2.0*ctx->H/ctx->omega_s*a + ctx->D;   J[1][0] = Pmax*PetscCosScalar(u[0]);

  ierr    = MatSetValues(B,2,rowcol,2,rowcol,&J[0][0],INSERT_VALUES);CHKERRQ(ierr);
  ierr    = VecRestoreArray(U,&u);CHKERRQ(ierr);
  ierr    = VecRestoreArray(Udot,&udot);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (A != B) {
    ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "RHSJacobianP"
static PetscErrorCode RHSJacobianP(TS ts,PetscReal t,Vec X,Mat A,void *ctx0)
{
  PetscErrorCode ierr;
  PetscInt       row[] = {0,1},col[]={0};
  PetscScalar    *x,J[2][1];
  PetscFunctionBeginUser;
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);

  J[0][0] = 0;
  J[1][0] = 1.;
  ierr    = MatSetValues(A,2,row,1,col,&J[0][0],INSERT_VALUES);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CostIntegrand"
static PetscErrorCode CostIntegrand(TS ts,PetscReal t,Vec U,Vec R,AppCtx *ctx)
{
  PetscErrorCode ierr;
  PetscScalar    *u,*r;
   
  PetscFunctionBegin; 
  ierr = VecGetArray(U,&u);CHKERRQ(ierr);
  ierr = VecGetArray(R,&r);CHKERRQ(ierr);
  
  r[0] = ctx->c*PetscPowScalarInt(PetscMax(0., u[0]-ctx->u_s),ctx->beta);CHKERRQ(ierr);
  /* r[0] = ctx->c*PetscPowScalarInt(u[0]-ctx->u_s,ctx->beta);CHKERRQ(ierr); */

  ierr = VecRestoreArray(R,&r);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DRDYFunction"
static PetscErrorCode DRDYFunction(TS ts,PetscReal t,Vec U,Vec *drdy,AppCtx *ctx)
{
  PetscErrorCode ierr;
  PetscScalar    *u,*ry;

  PetscFunctionBegin;
  ierr = VecGetArray(U,&u);CHKERRQ(ierr);
  ierr = VecGetArray(drdy[0],&ry);CHKERRQ(ierr);
  
  ry[0] = ctx->c*ctx->beta*PetscPowScalarInt(PetscMax(0., u[0]-ctx->u_s),ctx->beta-1);CHKERRQ(ierr);
  /* ry[0] = ctx->c*ctx->beta*PetscPowScalarInt(u[0]-ctx->u_s,ctx->beta-1);CHKERRQ(ierr); */
  ierr  = VecRestoreArray(drdy[0],&ry);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DRDPFunction"
static PetscErrorCode DRDPFunction(TS ts,PetscReal t,Vec U,Vec *drdp,AppCtx *ctx)
{
  PetscErrorCode ierr;
  PetscScalar    *u,*rp;

  PetscFunctionBegin;
  ierr = VecGetArray(U,&u);CHKERRQ(ierr);
  ierr = VecGetArray(drdp[0],&rp);CHKERRQ(ierr);

  rp[0] = 0.;
  ierr  = VecRestoreArray(drdp[0],&rp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeSensiP"
PetscErrorCode ComputeSensiP(Vec lambda,Vec lambdap,AppCtx *ctx)
{ 
  PetscErrorCode ierr;
  PetscScalar *x,*y,sensip;
  
  PetscFunctionBegin;
  ierr = VecGetArray(lambda,&x);CHKERRQ(ierr);
  ierr = VecGetArray(lambdap,&y);CHKERRQ(ierr);

  sensip = 1./PetscSqrtScalar(1.-(ctx->Pm/ctx->Pmax)*(ctx->Pm/ctx->Pmax))/ctx->Pmax*x[0]+y[0];  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n sensitivity wrt parameter pm: %.7f \n",sensip);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  TS             ts;            /* ODE integrator */
  Vec            U;             /* solution will be stored here */
  Mat            A;             /* Jacobian matrix */
  Mat            Jacp;          /* Jacobian matrix */
  PetscErrorCode ierr;
  PetscMPIInt    size;
  PetscInt       n = 2;
  AppCtx         ctx;
  PetscScalar    *u;
  PetscReal      du[2] = {0.0,0.0};
  PetscBool      ensemble = PETSC_FALSE,flg1,flg2;
  PetscReal      ftime;
  PetscInt       steps;
  PetscScalar    *x_ptr,*y_ptr;
  Vec            lambda[1],q,drdy[1],lambdap[1],drdp[1];
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Only for sequential runs");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Create necessary matrix and vectors
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,n,n,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = MatSetType(A,MATDENSE);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetVecs(A,&U,NULL);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&Jacp);CHKERRQ(ierr);
  ierr = MatSetSizes(Jacp,PETSC_DECIDE,PETSC_DECIDE,2,1);CHKERRQ(ierr);
  ierr = MatSetFromOptions(Jacp);CHKERRQ(ierr);
  ierr = MatSetUp(Jacp);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Set runtime options
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"Swing equation options","");CHKERRQ(ierr);
  {
    ctx.beta    = 2;
    ctx.c       = 10000.0;
    ctx.u_s     = 1.0;
    ctx.omega_s = 1.0;
    ctx.omega_b = 120.0*PETSC_PI;
    ctx.H       = 5.0;
    ierr        = PetscOptionsScalar("-Inertia","","",ctx.H,&ctx.H,NULL);CHKERRQ(ierr);
    ctx.D       = 5.0;
    ierr        = PetscOptionsScalar("-D","","",ctx.D,&ctx.D,NULL);CHKERRQ(ierr);
    ctx.E       = 1.1378;
    ctx.V       = 1.0;
    ctx.X       = 0.545;
    ctx.Pmax    = ctx.E*ctx.V/ctx.X;;
    ierr        = PetscOptionsScalar("-Pmax","","",ctx.Pmax,&ctx.Pmax,NULL);CHKERRQ(ierr);
    ctx.Pm      = 0.4;
    ierr        = PetscOptionsScalar("-Pm","","",ctx.Pm,&ctx.Pm,NULL);CHKERRQ(ierr);
    ctx.tf      = 0.1;
    ctx.tcl     = 0.2;
    ierr        = PetscOptionsReal("-tf","Time to start fault","",ctx.tf,&ctx.tf,NULL);CHKERRQ(ierr);
    ierr        = PetscOptionsReal("-tcl","Time to end fault","",ctx.tcl,&ctx.tcl,NULL);CHKERRQ(ierr);
    ierr        = PetscOptionsBool("-ensemble","Run ensemble of different initial conditions","",ensemble,&ensemble,NULL);CHKERRQ(ierr);
    if (ensemble) {
      ctx.tf      = -1;
      ctx.tcl     = -1;
    }

    ierr = VecGetArray(U,&u);CHKERRQ(ierr);
    u[0] = PetscAsinScalar(ctx.Pm/ctx.Pmax);
    u[1] = 1.0;
    ierr = PetscOptionsRealArray("-u","Initial solution","",u,&n,&flg1);CHKERRQ(ierr);
    n    = 2;
    ierr = PetscOptionsRealArray("-du","Perturbation in initial solution","",du,&n,&flg2);CHKERRQ(ierr);
    u[0] += du[0];
    u[1] += du[1];
    ierr = VecRestoreArray(U,&u);CHKERRQ(ierr);
    if (flg1 || flg2) {
      ctx.tf      = -1;
      ctx.tcl     = -1;
    }
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create timestepping solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetIFunction(ts,NULL,(TSIFunction) IFunction,&ctx);CHKERRQ(ierr);
  ierr = TSSetIJacobian(ts,A,A,(TSIJacobian)IJacobian,&ctx);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set initial conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSSetSolution(ts,U);CHKERRQ(ierr);

  /*
    Save trajectory of solution so that TSAdjointSolve() may be used
  */
  ierr = TSSetSaveTrajectory(ts);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set solver options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSSetDuration(ts,100000,10.0);CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,0.0,.01);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (ensemble) {
    for (du[1] = -2.5; du[1] <= .01; du[1] += .1) {
      ierr = VecGetArray(U,&u);CHKERRQ(ierr);
      u[0] = PetscAsinScalar(ctx.Pm/ctx.Pmax);
      u[1] = ctx.omega_s;
      u[0] += du[0];
      u[1] += du[1];
      ierr = VecRestoreArray(U,&u);CHKERRQ(ierr);
      ierr = TSSetInitialTimeStep(ts,0.0,.01);CHKERRQ(ierr);
      ierr = TSSolve(ts,U);CHKERRQ(ierr);
    }
  } else {
    ierr = TSSolve(ts,U);CHKERRQ(ierr);
  }

  ierr = TSGetSolveTime(ts,&ftime);CHKERRQ(ierr);
  ierr = TSGetTimeStepNumber(ts,&steps);CHKERRQ(ierr);
  ierr = VecView(U,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Adjoint model starts here
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MatGetVecs(A,&lambda[0],NULL);CHKERRQ(ierr);
  /*   Set initial conditions for the adjoint integration */
  ierr = VecGetArray(lambda[0],&y_ptr);CHKERRQ(ierr);
  y_ptr[0] = 0.0; y_ptr[1] = 0.0;
  ierr = VecRestoreArray(lambda[0],&y_ptr);CHKERRQ(ierr);
  ierr = TSAdjointSetSensitivity(ts,lambda,1);CHKERRQ(ierr);

  ierr = MatGetVecs(Jacp,&lambdap[0],NULL);CHKERRQ(ierr);
  ierr = VecGetArray(lambdap[0],&x_ptr);CHKERRQ(ierr);
  x_ptr[0] = -1.0;
  ierr = VecRestoreArray(lambdap[0],&x_ptr);CHKERRQ(ierr);
  ierr = TSAdjointSetSensitivityP(ts,lambdap,1);CHKERRQ(ierr);

  ierr = VecDuplicate(lambda[0],&drdy[0]);CHKERRQ(ierr);
  ierr = VecDuplicate(lambdap[0],&drdp[0]);CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_WORLD,1,&q);CHKERRQ(ierr);

  /*   Reset start time for the adjoint integration */
  ierr = TSSetTime(ts,ftime);CHKERRQ(ierr);

  /*   Set RHS Jacobian and number of steps for the adjoint integration */
  ierr = TSSetDuration(ts,steps,PETSC_DEFAULT);CHKERRQ(ierr);

  /*   Set RHS JacobianP */
  ierr = TSAdjointSetRHSJacobianP(ts,Jacp,RHSJacobianP,&ctx);CHKERRQ(ierr);

  ierr = TSAdjointSetCostIntegrand(ts,1,q,   (PetscErrorCode (*)(TS,PetscReal,Vec,Vec,void*))CostIntegrand,
                                        drdy,(PetscErrorCode (*)(TS,PetscReal,Vec,Vec*,void*))DRDYFunction,
                                        drdp,(PetscErrorCode (*)(TS,PetscReal,Vec,Vec*,void*))DRDPFunction,&ctx);CHKERRQ(ierr);

  ierr = TSAdjointSolve(ts);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n sensitivity wrt initial conditions: d[Psi(tf)]/d[phi0]  d[Psi(tf)]/d[omega0]\n");CHKERRQ(ierr);
  ierr = VecView(lambda[0],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecView(lambdap[0],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecView(q,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecGetArray(q,&x_ptr);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n cost function=%.15f\n",x_ptr[0]-ctx.Pm);CHKERRQ(ierr);
 
  ierr = ComputeSensiP(lambda[0],lambdap[0],&ctx);CHKERRQ(ierr);
/*
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n sensitivity wrt parameters: d[y(tf)]/d[mu]\n");CHKERRQ(ierr);
  ierr = VecView(lambdap[0],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
*/

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&Jacp);CHKERRQ(ierr);
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&lambda[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&drdy[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&drdp[0]);CHKERRQ(ierr);
  ierr = VecDestroy(&q);CHKERRQ(ierr);
  ierr = VecDestroy(&lambdap[0]);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return(0);
}
