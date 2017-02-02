/*
 * ylapack_cblas.h --
 *
 * Definitions for CBLAS functions.  This file is mostly a copy of CBLAS
 * standard header available at netlib (http://www.netlib.org/blas/) with
 * some additional macros for easier inclusion.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2011-2017 Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
 *
 * This file is part of YLAPACK.
 *
 * YLAPACK is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * YLAPACK is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with YLAPACK.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _CBLAS_H
#define _CBLAS_H 1
#include <stddef.h>

#undef _CBLAS_BEGIN_DECLS
#undef _CBLAS_END_DECLS
#ifdef __cplusplus
#define _CBLAS_BEGIN_DECLS extern "C" {
#define _CBLAS_END_DECLS }
#else
#define _CBLAS_BEGIN_DECLS           /* empty */
#define _CBLAS_END_DECLS             /* empty */
#endif

#ifndef CBLAS_INDEX
# define CBLAS_INDEX size_t  /* this may vary between platforms */
#endif
#ifndef CBLAS_INTEGER
# define CBLAS_INTEGER int
#endif

typedef enum CBLAS_ORDER     {CblasRowMajor    = 101,
                              CblasColMajor    = 102} CBLAS_ORDER;
typedef enum CBLAS_TRANSPOSE {CblasNoTrans     = 111,
                              CblasTrans       = 112,
                              CblasConjTrans   = 113,
                              CblasConjNoTrans = 114} CBLAS_TRANSPOSE;
typedef enum CBLAS_UPLO      {CblasUpper       = 121,
                              CblasLower       = 122} CBLAS_UPLO;
typedef enum CBLAS_DIAG      {CblasNonUnit     = 131,
                              CblasUnit        = 132} CBLAS_DIAG;
typedef enum CBLAS_SIDE      {CblasLeft        = 141,
                              CblasRight       = 142} CBLAS_SIDE;

_CBLAS_BEGIN_DECLS

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS functions (complex are recast as routines)
 * ===========================================================================
 */
float  cblas_sdsdot(const CBLAS_INTEGER N, const float alpha,
                    const float *X, const CBLAS_INTEGER incX,
                    const float *Y, const CBLAS_INTEGER incY);
double cblas_dsdot(const CBLAS_INTEGER N,
                   const float *X, const CBLAS_INTEGER incX,
                   const float *Y, const CBLAS_INTEGER incY);
float  cblas_sdot(const CBLAS_INTEGER N,
                  const float  *X, const CBLAS_INTEGER incX,
                  const float  *Y, const CBLAS_INTEGER incY);
double cblas_ddot(const CBLAS_INTEGER N,
                  const double *X, const CBLAS_INTEGER incX,
                  const double *Y, const CBLAS_INTEGER incY);

/*
 * Functions having prefixes Z and C only
 */
void   cblas_cdotu_sub(const CBLAS_INTEGER N,
                       const void *X, const CBLAS_INTEGER incX,
                       const void *Y, const CBLAS_INTEGER incY, void *dotu);
void   cblas_cdotc_sub(const CBLAS_INTEGER N,
                       const void *X, const CBLAS_INTEGER incX,
                       const void *Y, const CBLAS_INTEGER incY, void *dotc);

void   cblas_zdotu_sub(const CBLAS_INTEGER N,
                       const void *X, const CBLAS_INTEGER incX,
                       const void *Y, const CBLAS_INTEGER incY, void *dotu);
void   cblas_zdotc_sub(const CBLAS_INTEGER N,
                       const void *X, const CBLAS_INTEGER incX,
                       const void *Y, const CBLAS_INTEGER incY, void *dotc);


/*
 * Functions having prefixes S D SC DZ
 */
float  cblas_snrm2(const CBLAS_INTEGER N,
                   const float *X, const CBLAS_INTEGER incX);
float  cblas_sasum(const CBLAS_INTEGER N,
                   const float *X, const CBLAS_INTEGER incX);

double cblas_dnrm2(const CBLAS_INTEGER N,
                   const double *X, const CBLAS_INTEGER incX);
double cblas_dasum(const CBLAS_INTEGER N,
                   const double *X, const CBLAS_INTEGER incX);

float  cblas_scnrm2(const CBLAS_INTEGER N,
                    const void *X, const CBLAS_INTEGER incX);
float  cblas_scasum(const CBLAS_INTEGER N,
                    const void *X, const CBLAS_INTEGER incX);

double cblas_dznrm2(const CBLAS_INTEGER N,
                    const void *X, const CBLAS_INTEGER incX);
double cblas_dzasum(const CBLAS_INTEGER N,
                    const void *X, const CBLAS_INTEGER incX);


/*
 * Functions having standard 4 prefixes (S D C Z)
 */
CBLAS_INDEX cblas_isamax(const CBLAS_INTEGER N,
                         const float  *X, const CBLAS_INTEGER incX);
CBLAS_INDEX cblas_idamax(const CBLAS_INTEGER N,
                         const double *X, const CBLAS_INTEGER incX);
CBLAS_INDEX cblas_icamax(const CBLAS_INTEGER N,
                         const void   *X, const CBLAS_INTEGER incX);
CBLAS_INDEX cblas_izamax(const CBLAS_INTEGER N,
                         const void   *X, const CBLAS_INTEGER incX);

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS routines
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (s, d, c, z)
 */
void cblas_sswap(const CBLAS_INTEGER N,
                       float *X, const CBLAS_INTEGER incX,
                       float *Y, const CBLAS_INTEGER incY);
void cblas_scopy(const CBLAS_INTEGER N,
                 const float *X, const CBLAS_INTEGER incX,
                       float *Y, const CBLAS_INTEGER incY);
void cblas_saxpy(const CBLAS_INTEGER N, const float alpha,
                 const float *X, const CBLAS_INTEGER incX,
                       float *Y, const CBLAS_INTEGER incY);

void cblas_dswap(const CBLAS_INTEGER N,
                       double *X, const CBLAS_INTEGER incX,
                       double *Y, const CBLAS_INTEGER incY);
void cblas_dcopy(const CBLAS_INTEGER N,
                 const double *X, const CBLAS_INTEGER incX,
                       double *Y, const CBLAS_INTEGER incY);
void cblas_daxpy(const CBLAS_INTEGER N, const double alpha,
                 const double *X, const CBLAS_INTEGER incX,
                       double *Y, const CBLAS_INTEGER incY);

void cblas_cswap(const CBLAS_INTEGER N,
                       void *X, const CBLAS_INTEGER incX,
                       void *Y, const CBLAS_INTEGER incY);
void cblas_ccopy(const CBLAS_INTEGER N,
                 const void *X, const CBLAS_INTEGER incX,
                       void *Y, const CBLAS_INTEGER incY);
void cblas_caxpy(const CBLAS_INTEGER N, const void *alpha,
                 const void *X, const CBLAS_INTEGER incX,
                       void *Y, const CBLAS_INTEGER incY);

void cblas_zswap(const CBLAS_INTEGER N,
                       void *X, const CBLAS_INTEGER incX,
                       void *Y, const CBLAS_INTEGER incY);
void cblas_zcopy(const CBLAS_INTEGER N,
                 const void *X, const CBLAS_INTEGER incX,
                       void *Y, const CBLAS_INTEGER incY);
void cblas_zaxpy(const CBLAS_INTEGER N, const void *alpha,
                 const void *X, const CBLAS_INTEGER incX,
                       void *Y, const CBLAS_INTEGER incY);


/*
 * Routines with S and D prefix only
 */
void cblas_srotg(float *a, float *b, float *c, float *s);
void cblas_srotmg(float *d1, float *d2, float *b1, const float b2, float *P);
void cblas_srot(const CBLAS_INTEGER N,
                float *X, const CBLAS_INTEGER incX,
                float *Y, const CBLAS_INTEGER incY,
                const float c, const float s);
void cblas_srotm(const CBLAS_INTEGER N,
                 float *X, const CBLAS_INTEGER incX,
                 float *Y, const CBLAS_INTEGER incY, const float *P);

void cblas_drotg(double *a, double *b, double *c, double *s);
void cblas_drotmg(double *d1, double *d2, double *b1, const double b2,
                  double *P);
void cblas_drot(const CBLAS_INTEGER N,
                double *X, const CBLAS_INTEGER incX,
                double *Y, const CBLAS_INTEGER incY,
                const double c, const double  s);
void cblas_drotm(const CBLAS_INTEGER N,
                 double *X, const CBLAS_INTEGER incX,
                 double *Y, const CBLAS_INTEGER incY, const double *P);


/*
 * Routines with S D C Z CS and ZD prefixes
 */
void cblas_sscal(const CBLAS_INTEGER N, const float alpha,
                 float *X, const CBLAS_INTEGER incX);
void cblas_dscal(const CBLAS_INTEGER N, const double alpha,
                 double *X, const CBLAS_INTEGER incX);
void cblas_cscal(const CBLAS_INTEGER N, const void *alpha,
                 void *X, const CBLAS_INTEGER incX);
void cblas_zscal(const CBLAS_INTEGER N, const void *alpha,
                 void *X, const CBLAS_INTEGER incX);
void cblas_csscal(const CBLAS_INTEGER N, const float alpha,
                  void *X, const CBLAS_INTEGER incX);
void cblas_zdscal(const CBLAS_INTEGER N, const double alpha,
                  void *X, const CBLAS_INTEGER incX);

/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void cblas_sgemv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_INTEGER M, const CBLAS_INTEGER N,
                 const float alpha, const float *A, const CBLAS_INTEGER lda,
                 const float *X, const CBLAS_INTEGER incX, const float beta,
                 float *Y, const CBLAS_INTEGER incY);
void cblas_sgbmv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_INTEGER M, const CBLAS_INTEGER N,
                 const CBLAS_INTEGER KL, const CBLAS_INTEGER KU,
                 const float alpha, const float *A, const CBLAS_INTEGER lda,
                 const float *X, const CBLAS_INTEGER incX, const float beta,
                 float *Y, const CBLAS_INTEGER incY);
void cblas_strmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const float *A,
                 const CBLAS_INTEGER lda, float *X, const CBLAS_INTEGER incX);
void cblas_stbmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K, const float *A,
                 const CBLAS_INTEGER lda, float *X, const CBLAS_INTEGER incX);
void cblas_stpmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const float *Ap, float *X,
                 const CBLAS_INTEGER incX);
void cblas_strsv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const float *A,
                 const CBLAS_INTEGER lda, float *X, const CBLAS_INTEGER incX);
void cblas_stbsv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K, const float *A,
                 const CBLAS_INTEGER lda, float *X, const CBLAS_INTEGER incX);
void cblas_stpsv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const float *Ap, float *X,
                 const CBLAS_INTEGER incX);

void cblas_dgemv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_INTEGER M, const CBLAS_INTEGER N,
                 const double alpha, const double *A, const CBLAS_INTEGER lda,
                 const double *X, const CBLAS_INTEGER incX, const double beta,
                 double *Y, const CBLAS_INTEGER incY);
void cblas_dgbmv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_INTEGER M, const CBLAS_INTEGER N,
                 const CBLAS_INTEGER KL, const CBLAS_INTEGER KU,
                 const double alpha, const double *A, const CBLAS_INTEGER lda,
                 const double *X, const CBLAS_INTEGER incX, const double beta,
                 double *Y, const CBLAS_INTEGER incY);
void cblas_dtrmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const double *A,
                 const CBLAS_INTEGER lda, double *X, const CBLAS_INTEGER incX);
void cblas_dtbmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K,
                 const double *A, const CBLAS_INTEGER lda,
                 double *X, const CBLAS_INTEGER incX);
void cblas_dtpmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const double *Ap, double *X,
                 const CBLAS_INTEGER incX);
void cblas_dtrsv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const double *A,
                 const CBLAS_INTEGER lda, double *X, const CBLAS_INTEGER incX);
void cblas_dtbsv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K,
                 const double *A, const CBLAS_INTEGER lda,
                 double *X, const CBLAS_INTEGER incX);
void cblas_dtpsv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const double *Ap, double *X,
                 const CBLAS_INTEGER incX);

void cblas_cgemv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_INTEGER M, const CBLAS_INTEGER N,
                 const void *alpha, const void *A, const CBLAS_INTEGER lda,
                 const void *X, const CBLAS_INTEGER incX, const void *beta,
                 void *Y, const CBLAS_INTEGER incY);
void cblas_cgbmv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_INTEGER M, const CBLAS_INTEGER N,
                 const CBLAS_INTEGER KL, const CBLAS_INTEGER KU,
                 const void *alpha, const void *A, const CBLAS_INTEGER lda,
                 const void *X, const CBLAS_INTEGER incX, const void *beta,
                 void *Y, const CBLAS_INTEGER incY);
void cblas_ctrmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const void *A, const CBLAS_INTEGER lda,
                 void *X, const CBLAS_INTEGER incX);
void cblas_ctbmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K,
                 const void *A, const CBLAS_INTEGER lda,
                 void *X, const CBLAS_INTEGER incX);
void cblas_ctpmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const void *Ap, void *X,
                 const CBLAS_INTEGER incX);
void cblas_ctrsv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const void *A,
                 const CBLAS_INTEGER lda, void *X,
                 const CBLAS_INTEGER incX);
void cblas_ctbsv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K,
                 const void *A, const CBLAS_INTEGER lda,
                 void *X, const CBLAS_INTEGER incX);
void cblas_ctpsv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const void *Ap, void *X,
                 const CBLAS_INTEGER incX);

void cblas_zgemv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_INTEGER M, const CBLAS_INTEGER N,
                 const void *alpha, const void *A, const CBLAS_INTEGER lda,
                 const void *X, const CBLAS_INTEGER incX, const void *beta,
                 void *Y, const CBLAS_INTEGER incY);
void cblas_zgbmv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_INTEGER M, const CBLAS_INTEGER N,
                 const CBLAS_INTEGER KL, const CBLAS_INTEGER KU,
                 const void *alpha, const void *A, const CBLAS_INTEGER lda,
                 const void *X, const CBLAS_INTEGER incX, const void *beta,
                 void *Y, const CBLAS_INTEGER incY);
void cblas_ztrmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const void *A, const CBLAS_INTEGER lda,
                 void *X, const CBLAS_INTEGER incX);
void cblas_ztbmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K, const void *A,
                 const CBLAS_INTEGER lda,
                 void *X, const CBLAS_INTEGER incX);
void cblas_ztpmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const void *Ap, void *X,
                 const CBLAS_INTEGER incX);
void cblas_ztrsv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const void *A, const CBLAS_INTEGER lda,
                 void *X, const CBLAS_INTEGER incX);
void cblas_ztbsv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K, const void *A,
                 const CBLAS_INTEGER lda,
                 void *X, const CBLAS_INTEGER incX);
void cblas_ztpsv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const CBLAS_INTEGER N, const void *Ap, void *X,
                 const CBLAS_INTEGER incX);


/*
 * Routines with S and D prefixes only
 */
void cblas_ssymv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const float alpha, const float *A,
                 const CBLAS_INTEGER lda, const float *X,
                 const CBLAS_INTEGER incX, const float beta, float *Y,
                 const CBLAS_INTEGER incY);
void cblas_ssbmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K,
                 const float alpha, const float *A, const CBLAS_INTEGER lda,
                 const float *X, const CBLAS_INTEGER incX,
                 const float beta, float *Y, const CBLAS_INTEGER incY);
void cblas_sspmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const float alpha, const float *Ap,
                 const float *X, const CBLAS_INTEGER incX,
                 const float beta, float *Y, const CBLAS_INTEGER incY);
void cblas_sger(const CBLAS_ORDER order, const CBLAS_INTEGER M,
                const CBLAS_INTEGER N, const float alpha, const float *X,
                const CBLAS_INTEGER incX, const float *Y,
                const CBLAS_INTEGER incY, float *A, const CBLAS_INTEGER lda);
void cblas_ssyr(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                const CBLAS_INTEGER N, const float alpha, const float *X,
                const CBLAS_INTEGER incX, float *A, const CBLAS_INTEGER lda);
void cblas_sspr(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                const CBLAS_INTEGER N, const float alpha, const float *X,
                const CBLAS_INTEGER incX, float *Ap);
void cblas_ssyr2(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const float alpha, const float *X,
                 const CBLAS_INTEGER incX, const float *Y,
                 const CBLAS_INTEGER incY, float *A, const CBLAS_INTEGER lda);
void cblas_sspr2(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const float alpha, const float *X,
                 const CBLAS_INTEGER incX, const float *Y,
                 const CBLAS_INTEGER incY, float *A);

void cblas_dsymv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const double alpha, const double *A,
                 const CBLAS_INTEGER lda, const double *X,
                 const CBLAS_INTEGER incX,
                 const double beta, double *Y, const CBLAS_INTEGER incY);
void cblas_dsbmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K,
                 const double alpha, const double *A, const CBLAS_INTEGER lda,
                 const double *X, const CBLAS_INTEGER incX, const double beta,
                 double *Y, const CBLAS_INTEGER incY);
void cblas_dspmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const double alpha, const double *Ap,
                 const double *X, const CBLAS_INTEGER incX,
                 const double beta, double *Y, const CBLAS_INTEGER incY);
void cblas_dger(const CBLAS_ORDER order, const CBLAS_INTEGER M,
                const CBLAS_INTEGER N, const double alpha, const double *X,
                const CBLAS_INTEGER incX, const double *Y,
                const CBLAS_INTEGER incY, double *A, const CBLAS_INTEGER lda);
void cblas_dsyr(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                const CBLAS_INTEGER N, const double alpha, const double *X,
                const CBLAS_INTEGER incX, double *A, const CBLAS_INTEGER lda);
void cblas_dspr(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                const CBLAS_INTEGER N, const double alpha, const double *X,
                const CBLAS_INTEGER incX, double *Ap);
void cblas_dsyr2(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const double alpha, const double *X,
                 const CBLAS_INTEGER incX, const double *Y,
                 const CBLAS_INTEGER incY, double *A, const CBLAS_INTEGER lda);
void cblas_dspr2(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const double alpha, const double *X,
                 const CBLAS_INTEGER incX, const double *Y,
                 const CBLAS_INTEGER incY, double *A);


/*
 * Routines with C and Z prefixes only
 */
void cblas_chemv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const void *alpha, const void *A,
                 const CBLAS_INTEGER lda, const void *X,
                 const CBLAS_INTEGER incX, const void *beta, void *Y,
                 const CBLAS_INTEGER incY);
void cblas_chbmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K,
                 const void *alpha, const void *A, const CBLAS_INTEGER lda,
                 const void *X, const CBLAS_INTEGER incX,
                 const void *beta, void *Y, const CBLAS_INTEGER incY);
void cblas_chpmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const void *alpha, const void *Ap,
                 const void *X, const CBLAS_INTEGER incX,
                 const void *beta, void *Y, const CBLAS_INTEGER incY);
void cblas_cgeru(const CBLAS_ORDER order, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const void *alpha, const void *X,
                 const CBLAS_INTEGER incX, const void *Y,
                 const CBLAS_INTEGER incY, void *A, const CBLAS_INTEGER lda);
void cblas_cgerc(const CBLAS_ORDER order, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const void *alpha, const void *X,
                 const CBLAS_INTEGER incX, const void *Y,
                 const CBLAS_INTEGER incY, void *A, const CBLAS_INTEGER lda);
void cblas_cher(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                const CBLAS_INTEGER N, const float alpha, const void *X,
                const CBLAS_INTEGER incX, void *A, const CBLAS_INTEGER lda);
void cblas_chpr(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                const CBLAS_INTEGER N, const float alpha, const void *X,
                const CBLAS_INTEGER incX, void *A);
void cblas_cher2(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const void *alpha, const void *X,
                 const CBLAS_INTEGER incX, const void *Y,
                 const CBLAS_INTEGER incY, void *A, const CBLAS_INTEGER lda);
void cblas_chpr2(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const void *alpha, const void *X,
                 const CBLAS_INTEGER incX, const void *Y,
                 const CBLAS_INTEGER incY, void *Ap);

void cblas_zhemv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const void *alpha, const void *A,
                 const CBLAS_INTEGER lda, const void *X,
                 const CBLAS_INTEGER incX, const void *beta, void *Y,
                 const CBLAS_INTEGER incY);
void cblas_zhbmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K,
                 const void *alpha, const void *A, const CBLAS_INTEGER lda,
                 const void *X, const CBLAS_INTEGER incX, const void *beta,
                 void *Y, const CBLAS_INTEGER incY);
void cblas_zhpmv(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const void *alpha, const void *Ap,
                 const void *X, const CBLAS_INTEGER incX,
                 const void *beta, void *Y, const CBLAS_INTEGER incY);
void cblas_zgeru(const CBLAS_ORDER order, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const void *alpha, const void *X,
                 const CBLAS_INTEGER incX, const void *Y,
                 const CBLAS_INTEGER incY, void *A, const CBLAS_INTEGER lda);
void cblas_zgerc(const CBLAS_ORDER order, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const void *alpha, const void *X,
                 const CBLAS_INTEGER incX, const void *Y,
                 const CBLAS_INTEGER incY, void *A, const CBLAS_INTEGER lda);
void cblas_zher(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                const CBLAS_INTEGER N, const double alpha, const void *X,
                const CBLAS_INTEGER incX, void *A, const CBLAS_INTEGER lda);
void cblas_zhpr(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                const CBLAS_INTEGER N, const double alpha, const void *X,
                const CBLAS_INTEGER incX, void *A);
void cblas_zher2(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const void *alpha, const void *X,
                 const CBLAS_INTEGER incX, const void *Y,
                 const CBLAS_INTEGER incY, void *A, const CBLAS_INTEGER lda);
void cblas_zhpr2(const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
                 const CBLAS_INTEGER N, const void *alpha, const void *X,
                 const CBLAS_INTEGER incX, const void *Y,
                 const CBLAS_INTEGER incY, void *Ap);

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void cblas_sgemm(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_TRANSPOSE TransB, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K,
                 const float alpha, const float *A, const CBLAS_INTEGER lda,
                 const float *B, const CBLAS_INTEGER ldb, const float beta,
                 float *C, const CBLAS_INTEGER ldc);
void cblas_ssymm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const float alpha, const float *A,
                 const CBLAS_INTEGER lda, const float *B,
                 const CBLAS_INTEGER ldb, const float beta,
                 float *C, const CBLAS_INTEGER ldc);
void cblas_ssyrk(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE Trans, const CBLAS_INTEGER N,
                 const CBLAS_INTEGER K, const float alpha, const float *A,
                 const CBLAS_INTEGER lda, const float beta, float *C,
                 const CBLAS_INTEGER ldc);
void cblas_ssyr2k(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                  const CBLAS_TRANSPOSE Trans, const CBLAS_INTEGER N,
                  const CBLAS_INTEGER K, const float alpha, const float *A,
                  const CBLAS_INTEGER lda, const float *B,
                  const CBLAS_INTEGER ldb, const float beta, float *C,
                  const CBLAS_INTEGER ldc);
void cblas_strmm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_DIAG Diag, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const float alpha, const float *A,
                 const CBLAS_INTEGER lda, float *B, const CBLAS_INTEGER ldb);
void cblas_strsm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_DIAG Diag, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const float alpha, const float *A,
                 const CBLAS_INTEGER lda, float *B, const CBLAS_INTEGER ldb);

void cblas_dgemm(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_TRANSPOSE TransB, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K,
                 const double alpha, const double *A, const CBLAS_INTEGER lda,
                 const double *B, const CBLAS_INTEGER ldb, const double beta,
                 double *C, const CBLAS_INTEGER ldc);
void cblas_dsymm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const double alpha, const double *A,
                 const CBLAS_INTEGER lda, const double *B,
                 const CBLAS_INTEGER ldb, const double beta, double *C,
                 const CBLAS_INTEGER ldc);
void cblas_dsyrk(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE Trans, const CBLAS_INTEGER N,
                 const CBLAS_INTEGER K, const double alpha, const double *A,
                 const CBLAS_INTEGER lda, const double beta, double *C,
                 const CBLAS_INTEGER ldc);
void cblas_dsyr2k(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                  const CBLAS_TRANSPOSE Trans, const CBLAS_INTEGER N,
                  const CBLAS_INTEGER K, const double alpha, const double *A,
                  const CBLAS_INTEGER lda, const double *B,
                  const CBLAS_INTEGER ldb, const double beta, double *C,
                  const CBLAS_INTEGER ldc);
void cblas_dtrmm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_DIAG Diag, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const double alpha, const double *A,
                 const CBLAS_INTEGER lda, double *B, const CBLAS_INTEGER ldb);
void cblas_dtrsm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_DIAG Diag, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const double alpha, const double *A,
                 const CBLAS_INTEGER lda, double *B, const CBLAS_INTEGER ldb);

void cblas_cgemm(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_TRANSPOSE TransB, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K,
                 const void *alpha, const void *A, const CBLAS_INTEGER lda,
                 const void *B, const CBLAS_INTEGER ldb, const void *beta,
                 void *C, const CBLAS_INTEGER ldc);
void cblas_csymm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const void *alpha, const void *A,
                 const CBLAS_INTEGER lda, const void *B,
                 const CBLAS_INTEGER ldb, const void *beta, void *C,
                 const CBLAS_INTEGER ldc);
void cblas_csyrk(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE Trans, const CBLAS_INTEGER N,
                 const CBLAS_INTEGER K, const void *alpha, const void *A,
                 const CBLAS_INTEGER lda, const void *beta, void *C,
                 const CBLAS_INTEGER ldc);
void cblas_csyr2k(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                  const CBLAS_TRANSPOSE Trans, const CBLAS_INTEGER N,
                  const CBLAS_INTEGER K, const void *alpha, const void *A,
                  const CBLAS_INTEGER lda, const void *B,
                  const CBLAS_INTEGER ldb, const void *beta, void *C,
                  const CBLAS_INTEGER ldc);
void cblas_ctrmm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_DIAG Diag, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const void *alpha, const void *A,
                 const CBLAS_INTEGER lda, void *B, const CBLAS_INTEGER ldb);
void cblas_ctrsm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_DIAG Diag, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const void *alpha, const void *A,
                 const CBLAS_INTEGER lda, void *B, const CBLAS_INTEGER ldb);

void cblas_zgemm(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_TRANSPOSE TransB, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const CBLAS_INTEGER K,
                 const void *alpha, const void *A, const CBLAS_INTEGER lda,
                 const void *B, const CBLAS_INTEGER ldb, const void *beta,
                 void *C, const CBLAS_INTEGER ldc);
void cblas_zsymm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N,
                 const void *alpha, const void *A, const CBLAS_INTEGER lda,
                 const void *B, const CBLAS_INTEGER ldb, const void *beta,
                 void *C, const CBLAS_INTEGER ldc);
void cblas_zsyrk(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE Trans, const CBLAS_INTEGER N,
                 const CBLAS_INTEGER K,
                 const void *alpha, const void *A, const CBLAS_INTEGER lda,
                 const void *beta, void *C, const CBLAS_INTEGER ldc);
void cblas_zsyr2k(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                  const CBLAS_TRANSPOSE Trans, const CBLAS_INTEGER N,
                  const CBLAS_INTEGER K,
                  const void *alpha, const void *A, const CBLAS_INTEGER lda,
                  const void *B, const CBLAS_INTEGER ldb, const void *beta,
                  void *C, const CBLAS_INTEGER ldc);
void cblas_ztrmm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_DIAG Diag, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N,
                 const void *alpha, const void *A, const CBLAS_INTEGER lda,
                 void *B, const CBLAS_INTEGER ldb);
void cblas_ztrsm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_DIAG Diag, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N,
                 const void *alpha, const void *A, const CBLAS_INTEGER lda,
                 void *B, const CBLAS_INTEGER ldb);


/*
 * Routines with prefixes C and Z only
 */
void cblas_chemm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const void *alpha, const void *A,
                 const CBLAS_INTEGER lda, const void *B,
                 const CBLAS_INTEGER ldb, const void *beta, void *C,
                 const CBLAS_INTEGER ldc);
void cblas_cherk(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE Trans, const CBLAS_INTEGER N,
                 const CBLAS_INTEGER K, const float alpha, const void *A,
                 const CBLAS_INTEGER lda, const float beta, void *C,
                 const CBLAS_INTEGER ldc);
void cblas_cher2k(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                  const CBLAS_TRANSPOSE Trans, const CBLAS_INTEGER N,
                  const CBLAS_INTEGER K, const void *alpha, const void *A,
                  const CBLAS_INTEGER lda, const void *B,
                  const CBLAS_INTEGER ldb, const float beta, void *C,
                  const CBLAS_INTEGER ldc);
void cblas_zhemm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_INTEGER M,
                 const CBLAS_INTEGER N, const void *alpha, const void *A,
                 const CBLAS_INTEGER lda, const void *B,
                 const CBLAS_INTEGER ldb, const void *beta, void *C,
                 const CBLAS_INTEGER ldc);
void cblas_zherk(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE Trans, const CBLAS_INTEGER N,
                 const CBLAS_INTEGER K, const double alpha, const void *A,
                 const CBLAS_INTEGER lda, const double beta, void *C,
                 const CBLAS_INTEGER ldc);
void cblas_zher2k(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                  const CBLAS_TRANSPOSE Trans, const CBLAS_INTEGER N,
                  const CBLAS_INTEGER K, const void *alpha, const void *A,
                  const CBLAS_INTEGER lda, const void *B,
                  const CBLAS_INTEGER ldb, const double beta, void *C,
                  const CBLAS_INTEGER ldc);


void cblas_xerbla(CBLAS_INTEGER p, const char *rout, const char *form, ...);

_CBLAS_END_DECLS

#endif /* _CBLAS_H */
