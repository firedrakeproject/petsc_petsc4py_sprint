#ifndef lint
static char vcid[] = "$Id: bdfact.c,v 1.38 1996/08/15 12:47:45 bsmith Exp bsmith $";
#endif

/* Block diagonal matrix format - factorization and triangular solves */

#include "src/mat/impls/bdiag/seq/bdiag.h"
#include "src/vec/vecimpl.h"
#include "src/inline/ilu.h"

int MatILUFactorSymbolic_SeqBDiag(Mat A,IS isrow,IS iscol,double f,
                                  int levels,Mat *B)
{
  Mat_SeqBDiag *a = (Mat_SeqBDiag *) A->data;
  PetscTruth   idn;
  int          ierr;

  if (a->m != a->n) SETERRQ(1,"MatILUFactorSymbolic_SeqBDiag:Matrix must be square");
  if (isrow) {
    ierr = ISIdentity(isrow,&idn); CHKERRQ(ierr);
    if (!idn) SETERRQ(1,"MatILUFactorSymbolic_SeqBDiag:Only identity row permutation supported");
  }
  if (iscol) {
    ierr = ISIdentity(iscol,&idn); CHKERRQ(ierr);
    if (!idn) SETERRQ(1,"MatILUFactorSymbolic_SeqBDiag:Only identity column permutation supported");
  }
  if (levels != 0)
    SETERRQ(1,"MatLUFactorSymbolic_SeqBDiag:Only ILU(0) is supported");
  ierr = MatConvert(A,MATSAME,B); CHKERRQ(ierr);

  /* Must set to zero for repeated calls with different nonzero structure */
  (*B)->factor = 0;
  return 0;
}

int MatILUFactor_SeqBDiag(Mat A,IS isrow,IS iscol,double f,int level)
{
  Mat_SeqBDiag *a = (Mat_SeqBDiag *) A->data;
  PetscTruth   idn;
  int          ierr;

  /* For now, no fill is allocated in symbolic factorization phase, so we
     directly use the input matrix for numeric factorization. */
  if (a->m != a->n) SETERRQ(1,"MatILUFactor_SeqBDiag:Matrix must be square");
  if (isrow) {
    ierr = ISIdentity(isrow,&idn); CHKERRQ(ierr);
    if (!idn) SETERRQ(1,"MatILUFactor_SeqBDiag:Only identity row permutation supported");
  }
  if (iscol) {
    ierr = ISIdentity(iscol,&idn); CHKERRQ(ierr);
    if (!idn) SETERRQ(1,"MatILUFactor_SeqBDiag:Only identity column permutation supported");
  }
  if (level != 0) SETERRQ(1,"MatILUFactor_SeqBDiag:Only ILU(0) is supported");
  return MatLUFactorNumeric(A,&A);
}

/* --------------------------------------------------------------------------*/
int MatLUFactorNumeric_SeqBDiag_N(Mat A,Mat *B)
{
  Mat          C = *B;
  Mat_SeqBDiag *a = (Mat_SeqBDiag *) C->data, *a1 = (Mat_SeqBDiag *) A->data;
  int          k, d, d2, dgk, elim_row, elim_col, bs = a->bs, knb, knb2, bs2 = bs*bs;
  int          dnum,nd = a->nd, mblock = a->mblock, nblock = a->nblock, ierr;
  int          *diag = a->diag,  m = a->m, mainbd = a->mainbd, *dgptr, len, i;
  Scalar       **dv = a->diagv, *dd = dv[mainbd], *v_work;
  Scalar       *multiplier;

  /* Copy input matrix to factored matrix if we've already factored the
     matrix before AND the nonzero structure remains the same.  This is done
     in symbolic factorization the first time through, but there's no symbolic
     factorization for successive calls with same matrix sparsity structure. */
  if (C->factor == FACTOR_LU) {
    for (i=0; i<a->nd; i++) {
      len = a->bdlen[i] * bs2 * sizeof(Scalar);
      d   = diag[i];
      if (d > 0) PetscMemcpy(dv[i]+bs2*d,a1->diagv[i]+bs2*d,len);
      else       PetscMemcpy(dv[i],a1->diagv[i],len);
    }
  }

  if (!a->pivot) {
    a->pivot = (int *) PetscMalloc(m*sizeof(int)); CHKPTRQ(a->pivot);
    PLogObjectMemory(C,m*sizeof(int));
  }
  v_work = (Scalar *) PetscMalloc((bs2+bs)*sizeof(Scalar));CHKPTRQ(v_work);
  multiplier = v_work + bs;
  dgptr = (int *) PetscMalloc((mblock+nblock)*sizeof(int)); CHKPTRQ(dgptr);
  PetscMemzero(dgptr,(mblock+nblock)*sizeof(int));
  for ( k=0; k<nd; k++ ) dgptr[diag[k]+mblock] = k+1;
  for ( k=0; k<mblock; k++ ) { /* k = block pivot_row */
    knb = k*bs; knb2 = knb*bs;
    /* invert the diagonal block */
    Kernel_A_gets_inverse_A(bs,dd+knb2,a->pivot+knb,v_work);
    for ( d=mainbd-1; d>=0; d-- ) {
      elim_row = k + diag[d];
      if (elim_row < mblock) { /* sweep down */
        /* dv[d][knb2]: test if entire block is zero? */
        Kernel_A_gets_A_times_B(bs,&dv[d][elim_row*bs2],dd+knb2,multiplier); 
        for ( d2=d+1; d2<nd; d2++ ) {
          elim_col = elim_row - diag[d2];
          if (elim_col >=0 && elim_col < nblock) {
            dgk = k - elim_col;
            if ((dnum = dgptr[dgk+mblock])) {
              Kernel_A_gets_A_minus_B_times_C(bs,&dv[d2][elim_row*bs2],
                             &dv[d][elim_row*bs2],&dv[dnum-1][knb2]);
            }
          }
        }
      }
    }
  }
  PetscFree(dgptr); PetscFree(v_work);
  C->factor = FACTOR_LU;
  return 0;
}

int MatLUFactorNumeric_SeqBDiag_1(Mat A,Mat *B)
{
  Mat          C = *B;
  Mat_SeqBDiag *a = (Mat_SeqBDiag *) C->data, *a1 = (Mat_SeqBDiag *) A->data;
  int          k, d, d2, dgk, elim_row, elim_col, dnum, nd = a->nd, i, len;
  int          *diag = a->diag, n = a->n, m = a->m, mainbd = a->mainbd, *dgptr;
  Scalar       **dv = a->diagv, *dd = dv[mainbd], mult;

  /* Copy input matrix to factored matrix if we've already factored the
     matrix before AND the nonzero structure remains the same.  This is done
     in symbolic factorization the first time through, but there's no symbolic
     factorization for successive calls with same matrix sparsity structure. */
  if (C->factor == FACTOR_LU) {
    for (i=0; i<nd; i++) {
      len = a->bdlen[i] * sizeof(Scalar);
      d   = diag[i];
      if (d > 0) PetscMemcpy(dv[i]+d,a1->diagv[i]+d,len);
      else       PetscMemcpy(dv[i],a1->diagv[i],len);
    }
  }

  dgptr = (int *) PetscMalloc((m+n)*sizeof(int)); CHKPTRQ(dgptr);
  PetscMemzero(dgptr,(m+n)*sizeof(int));
  for ( k=0; k<nd; k++ ) dgptr[diag[k]+m] = k+1;
  for ( k=0; k<m; k++ ) { /* k = pivot_row */
    dd[k] = 1.0/dd[k];
    for ( d=mainbd-1; d>=0; d-- ) {
      elim_row = k + diag[d];
      if (elim_row < m) { /* sweep down */
        if (dv[d][elim_row] != 0) {
          dv[d][elim_row] *= dd[k];
          mult = dv[d][elim_row];
          for ( d2=d+1; d2<nd; d2++ ) {
            elim_col = elim_row - diag[d2];
            dgk = k - elim_col;
            if (elim_col >=0 && elim_col < n) {
              if ((dnum = dgptr[dgk+m])) {
                dv[d2][elim_row] -= mult * dv[dnum-1][k];
              }
            }
          }
        }
      }
    }
  }
  PetscFree(dgptr);
  C->factor = FACTOR_LU;
  return 0;
}

/* -----------------------------------------------------------------*/

int MatSolve_SeqBDiag_1(Mat A,Vec xx,Vec yy)
{
  Mat_SeqBDiag *a = (Mat_SeqBDiag *) A->data;
  int          i, d, loc, mainbd = a->mainbd;
  int          n = a->n, m = a->m, *diag = a->diag, col;
  Scalar       *x, *y, *dd = a->diagv[mainbd], sum, **dv = a->diagv;

  VecGetArray_Fast(xx,x);
  VecGetArray_Fast(yy,y);
  /* forward solve the lower triangular part */
  for (i=0; i<m; i++) {
    sum = x[i];
    for (d=0; d<mainbd; d++) {
      loc = i - diag[d];
      if (loc >= 0) sum -= dv[d][i] * y[loc];
    }
    y[i] = sum;
  }
  /* backward solve the upper triangular part */
  for ( i=m-1; i>=0; i-- ) {
    sum = y[i];
    for (d=mainbd+1; d<a->nd; d++) {
      col = i - diag[d];
      if (col < n) sum -= dv[d][i] * y[col];
    }
    y[i] = sum*dd[i];
  }
  PLogFlops(2*a->nz - a->n);
  return 0;
}

int MatSolve_SeqBDiag_2(Mat A,Vec xx,Vec yy)
{
  Mat_SeqBDiag *a = (Mat_SeqBDiag *) A->data;
  int          i, d, loc, mainbd = a->mainbd;
  int          mblock = a->mblock, nblock = a->nblock, inb, inb2;
  int          m = a->m, *diag = a->diag, col;
  Scalar       *x, *y, *dd = a->diagv[mainbd], **dv = a->diagv,*dvt;
  Scalar       w0,w1,sum0,sum1;

  VecGetArray_Fast(xx,x);
  VecGetArray_Fast(yy,y);
  PetscMemcpy(y,x,m*sizeof(Scalar));

  /* forward solve the lower triangular part */
  if (mainbd != 0) {
    inb = 0;
    for (i=0; i<mblock; i++) {
      sum0 = sum1 = 0.0;
      for (d=0; d<mainbd; d++) {
        loc = 2*(i - diag[d]);
        if (loc >= 0) {
          dvt = &dv[d][4*i]; 
          w0 = y[loc]; w1 = y[loc+1];
          sum0 += dvt[0]*w0 + dvt[2]*w1;
          sum1 += dvt[1]*w0 + dvt[3]*w1;
        }
      }
      y[inb] -= sum0; y[inb+1] -= sum1; 

      inb += 2;
    }
  }
  /* backward solve the upper triangular part */
  inb = 2*(mblock-1); inb2 = 2*inb;
  for ( i=mblock-1; i>=0; i-- ) {
    sum0 = y[inb]; sum1 = y[inb+1];
    for (d=mainbd+1; d<a->nd; d++) {
      col = 2*(i - diag[d]);
      if (col < 2*nblock) {
        dvt = &dv[d][4*i]; 
        w0 = y[col]; w1 = y[col+1];
        sum0 -= dvt[0]*w0 + dvt[2]*w1;
        sum1 -= dvt[1]*w0 + dvt[3]*w1;
      }
    }
    dvt = dd+inb2;
    y[inb]   = dvt[0]*sum0 + dvt[2]*sum1;
    y[inb+1] = dvt[1]*sum0 + dvt[3]*sum1;
    inb -= 2; inb2 -= 4;
  }
  PLogFlops(2*a->nz - a->n);
  return 0;
}

int MatSolve_SeqBDiag_3(Mat A,Vec xx,Vec yy)
{
  Mat_SeqBDiag *a = (Mat_SeqBDiag *) A->data;
  int          i, d, loc, mainbd = a->mainbd;
  int          mblock = a->mblock, nblock = a->nblock, inb, inb2;
  int          m = a->m, *diag = a->diag, col;
  Scalar       *x, *y, *dd = a->diagv[mainbd], **dv = a->diagv,*dvt;
  Scalar       w0,w1,w2,sum0,sum1,sum2;

  VecGetArray_Fast(xx,x);
  VecGetArray_Fast(yy,y);
  PetscMemcpy(y,x,m*sizeof(Scalar));

  /* forward solve the lower triangular part */
  if (mainbd != 0) {
    inb = 0;
    for (i=0; i<mblock; i++) {
      sum0 = sum1 = sum2 = 0.0;
      for (d=0; d<mainbd; d++) {
        loc = 3*(i - diag[d]);
        if (loc >= 0) {
          dvt = &dv[d][9*i]; 
          w0 = y[loc]; w1 = y[loc+1]; w2 = y[loc+2];
          sum0 += dvt[0]*w0 + dvt[3]*w1 + dvt[6]*w2;
          sum1 += dvt[1]*w0 + dvt[4]*w1 + dvt[7]*w2;
          sum2 += dvt[2]*w0 + dvt[5]*w1 + dvt[8]*w2;
        }
      }
      y[inb] -= sum0; y[inb+1] -= sum1; y[inb+2] -= sum2;
      inb += 3;
    }
  }
  /* backward solve the upper triangular part */
  inb = 3*(mblock-1); inb2 = 3*inb;
  for ( i=mblock-1; i>=0; i-- ) {
    sum0 = y[inb]; sum1 = y[inb+1]; sum2 =  y[inb+2];
    for (d=mainbd+1; d<a->nd; d++) {
      col = 3*(i - diag[d]);
      if (col < 3*nblock) {
        dvt = &dv[d][9*i]; 
        w0 = y[col]; w1 = y[col+1];w2 = y[col+2];
        sum0 -= dvt[0]*w0 + dvt[3]*w1 + dvt[6]*w2;
        sum1 -= dvt[1]*w0 + dvt[4]*w1 + dvt[7]*w2;
        sum2 -= dvt[2]*w0 + dvt[5]*w1 + dvt[8]*w2;
      }
    }
    dvt = dd+inb2;
    y[inb]   = dvt[0]*sum0 + dvt[3]*sum1 + dvt[6]*sum2;
    y[inb+1] = dvt[1]*sum0 + dvt[4]*sum1 + dvt[7]*sum2;
    y[inb+2] = dvt[2]*sum0 + dvt[5]*sum1 + dvt[8]*sum2;
    inb -= 3; inb2 -= 9;
  }
  PLogFlops(2*a->nz - a->n);
  return 0;
}

int MatSolve_SeqBDiag_4(Mat A,Vec xx,Vec yy)
{
  Mat_SeqBDiag *a = (Mat_SeqBDiag *) A->data;
  int          i, d, loc, mainbd = a->mainbd;
  int          mblock = a->mblock, nblock = a->nblock, inb, inb2;
  int          m = a->m, *diag = a->diag, col;
  Scalar       *x, *y, *dd = a->diagv[mainbd], **dv = a->diagv,*dvt;
  Scalar       w0,w1,w2,w3,sum0,sum1,sum2,sum3;

  VecGetArray_Fast(xx,x);
  VecGetArray_Fast(yy,y);
  PetscMemcpy(y,x,m*sizeof(Scalar));

  /* forward solve the lower triangular part */
  if (mainbd != 0) {
    inb = 0;
    for (i=0; i<mblock; i++) {
      sum0 = sum1 = sum2 = sum3 = 0.0;
      for (d=0; d<mainbd; d++) {
        loc = 4*(i - diag[d]);
        if (loc >= 0) {
          dvt = &dv[d][16*i]; 
          w0 = y[loc]; w1 = y[loc+1]; w2 = y[loc+2];w3 = y[loc+3];
          sum0 += dvt[0]*w0 + dvt[4]*w1 + dvt[8]*w2  + dvt[12]*w3;
          sum1 += dvt[1]*w0 + dvt[5]*w1 + dvt[9]*w2  + dvt[13]*w3;
          sum2 += dvt[2]*w0 + dvt[6]*w1 + dvt[10]*w2 + dvt[14]*w3;
          sum3 += dvt[3]*w0 + dvt[7]*w1 + dvt[11]*w2 + dvt[15]*w3;
        }
      }
      y[inb] -= sum0; y[inb+1] -= sum1; y[inb+2] -= sum2;y[inb+3] -= sum3; 
      inb += 4;
    }
  }
  /* backward solve the upper triangular part */
  inb = 4*(mblock-1); inb2 = 4*inb;
  for ( i=mblock-1; i>=0; i-- ) {
    sum0 = y[inb]; sum1 = y[inb+1]; sum2 =  y[inb+2]; sum3 =  y[inb+3];
    for (d=mainbd+1; d<a->nd; d++) {
      col = 4*(i - diag[d]);
      if (col < 4*nblock) {
        dvt = &dv[d][16*i]; 
        w0 = y[col]; w1 = y[col+1];w2 = y[col+2];w3 = y[col+3];
        sum0 -= dvt[0]*w0 + dvt[4]*w1 + dvt[8]*w2  + dvt[12]*w3;
        sum1 -= dvt[1]*w0 + dvt[5]*w1 + dvt[9]*w2  + dvt[13]*w3;
        sum2 -= dvt[2]*w0 + dvt[6]*w1 + dvt[10]*w2 + dvt[14]*w3;
        sum3 -= dvt[3]*w0 + dvt[7]*w1 + dvt[11]*w2 + dvt[15]*w3;
      }
    }
    dvt = dd+inb2;
    y[inb]   = dvt[0]*sum0 + dvt[4]*sum1 + dvt[8]*sum2 + dvt[12]*sum3;
    y[inb+1] = dvt[1]*sum0 + dvt[5]*sum1 + dvt[9]*sum2 + dvt[13]*sum3;
    y[inb+2] = dvt[2]*sum0 + dvt[6]*sum1 + dvt[10]*sum2 + dvt[14]*sum3;
    y[inb+3] = dvt[3]*sum0 + dvt[7]*sum1 + dvt[11]*sum2 + dvt[15]*sum3;
    inb -= 4; inb2 -= 16;
  }
  PLogFlops(2*a->nz - a->n);
  return 0;
}

int MatSolve_SeqBDiag_5(Mat A,Vec xx,Vec yy)
{
  Mat_SeqBDiag *a = (Mat_SeqBDiag *) A->data;
  int          i, d, loc, mainbd = a->mainbd;
  int          mblock = a->mblock, nblock = a->nblock, inb, inb2;
  int          m = a->m, *diag = a->diag, col;
  Scalar       *x, *y, *dd = a->diagv[mainbd], **dv = a->diagv,*dvt;
  Scalar       w0,w1,w2,w3,w4,sum0,sum1,sum2,sum3,sum4;

  VecGetArray_Fast(xx,x);
  VecGetArray_Fast(yy,y);
  PetscMemcpy(y,x,m*sizeof(Scalar));

  /* forward solve the lower triangular part */
  if (mainbd != 0) {
    inb = 0;
    for (i=0; i<mblock; i++) {
      sum0 = sum1 = sum2 = sum3 = sum4 = 0.0;
      for (d=0; d<mainbd; d++) {
        loc = 5*(i - diag[d]);
        if (loc >= 0) {
          dvt = &dv[d][25*i]; 
          w0 = y[loc]; w1 = y[loc+1]; w2 = y[loc+2];w3 = y[loc+3];w4 = y[loc+4];
          sum0 += dvt[0]*w0 + dvt[5]*w1 + dvt[10]*w2 + dvt[15]*w3 + dvt[20]*w4;
          sum1 += dvt[1]*w0 + dvt[6]*w1 + dvt[11]*w2 + dvt[16]*w3 + dvt[21]*w4;
          sum2 += dvt[2]*w0 + dvt[7]*w1 + dvt[12]*w2 + dvt[17]*w3 + dvt[22]*w4;
          sum3 += dvt[3]*w0 + dvt[8]*w1 + dvt[13]*w2 + dvt[18]*w3 + dvt[23]*w4;
          sum4 += dvt[4]*w0 + dvt[9]*w1 + dvt[14]*w2 + dvt[19]*w3 + dvt[24]*w4;
        }
      }
      y[inb]   -= sum0; y[inb+1] -= sum1; y[inb+2] -= sum2;y[inb+3] -= sum3; 
      y[inb+4] -= sum4;
      inb += 5;
    }
  }
  /* backward solve the upper triangular part */
  inb = 5*(mblock-1); inb2 = 5*inb;
  for ( i=mblock-1; i>=0; i-- ) {
    sum0 = y[inb];sum1 = y[inb+1];sum2 = y[inb+2];sum3 = y[inb+3];sum4 = y[inb+4];
    for (d=mainbd+1; d<a->nd; d++) {
      col = 5*(i - diag[d]);
      if (col < 5*nblock) {
        dvt = &dv[d][25*i]; 
        w0 = y[col]; w1 = y[col+1];w2 = y[col+2];w3 = y[col+3];w4 = y[col+4];
        sum0 -= dvt[0]*w0 + dvt[5]*w1 + dvt[10]*w2 + dvt[15]*w3 + dvt[20]*w4;
        sum1 -= dvt[1]*w0 + dvt[6]*w1 + dvt[11]*w2 + dvt[16]*w3 + dvt[21]*w4;
        sum2 -= dvt[2]*w0 + dvt[7]*w1 + dvt[12]*w2 + dvt[17]*w3 + dvt[22]*w4;
        sum3 -= dvt[3]*w0 + dvt[8]*w1 + dvt[13]*w2 + dvt[18]*w3 + dvt[23]*w4;
        sum4 -= dvt[4]*w0 + dvt[9]*w1 + dvt[14]*w2 + dvt[19]*w3 + dvt[24]*w4;
      }
    }
    dvt = dd+inb2;
    y[inb]   = dvt[0]*sum0 + dvt[5]*sum1 + dvt[10]*sum2 + dvt[15]*sum3 
               + dvt[20]*sum4;
    y[inb+1] = dvt[1]*sum0 + dvt[6]*sum1 + dvt[11]*sum2 + dvt[16]*sum3 
               + dvt[21]*sum4;
    y[inb+2] = dvt[2]*sum0 + dvt[7]*sum1 + dvt[12]*sum2 + dvt[17]*sum3 
               + dvt[22]*sum4;
    y[inb+3] = dvt[3]*sum0 + dvt[8]*sum1 + dvt[13]*sum2 + dvt[18]*sum3 
               + dvt[23]*sum4;
    y[inb+4] = dvt[4]*sum0 + dvt[9]*sum1 + dvt[14]*sum2 + dvt[19]*sum3 
               + dvt[24]*sum4;
    inb -= 5; inb2 -= 25;
  }
  PLogFlops(2*a->nz - a->n);
  return 0;
}

int MatSolve_SeqBDiag_N(Mat A,Vec xx,Vec yy)
{
  Mat_SeqBDiag *a = (Mat_SeqBDiag *) A->data;
  int          i, d, loc, mainbd = a->mainbd;
  int          mblock = a->mblock, nblock = a->nblock, inb, inb2;
  int          bs = a->bs, m = a->m, *diag = a->diag, col, bs2 = bs*bs;
  Scalar       *x, *y, *dd = a->diagv[mainbd], **dv = a->diagv;
  Scalar       work[25];

  VecGetArray_Fast(xx,x);
  VecGetArray_Fast(yy,y);
  if (bs > 25) SETERRQ(1,"Blocks must be smaller then 25");
  PetscMemcpy(y,x,m*sizeof(Scalar));

  /* forward solve the lower triangular part */
  if (mainbd != 0) {
    inb = 0;
    for (i=0; i<mblock; i++) {
      for (d=0; d<mainbd; d++) {
        loc = i - diag[d];
        if (loc >= 0) {
          Kernel_v_gets_v_minus_A_times_w(bs,y+inb,&dv[d][i*bs2],y+loc*bs);
        }
      }
      inb += bs;
    }
  }
  /* backward solve the upper triangular part */
  inb = bs*(mblock-1); inb2 = inb*bs;
  for ( i=mblock-1; i>=0; i-- ) {
    for (d=mainbd+1; d<a->nd; d++) {
      col = i - diag[d];
      if (col < nblock) {
        Kernel_v_gets_v_minus_A_times_w(bs,y+inb,&dv[d][inb2],y+col*bs);
      }
    }
    Kernel_w_gets_A_times_v(bs,y+inb,dd+inb2,work);  
    PetscMemcpy(y+inb,work,bs*sizeof(Scalar));
    inb -= bs; inb2 -= bs2;
  }
  PLogFlops(2*a->nz - a->n);
  return 0;
}






