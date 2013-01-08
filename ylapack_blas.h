/* 
 * ylapack_blas.h --
 *
 * Definitions for BLAS functions.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2011 Éric Thiébaut <thiebaut@obs.univ-lyon1.fr>
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can use, modify
 * and/or redistribute the software under the terms of the CeCILL-C license as
 * circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited liability.
 *
 * In this respect, the user's attention is drawn to the risks associated with
 * loading, using, modifying and/or developing or reproducing the software by
 * the user in light of its specific status of free software, that may mean
 * that it is complicated to manipulate, and that also therefore means that it
 * is reserved for developers and experienced professionals having in-depth
 * computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling
 * the security of their systems and/or data to be ensured and, more
 * generally, to use and operate it in the same conditions as regards
 * security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 *
 *-----------------------------------------------------------------------------
 */

#ifndef _YLAPACK_H
# error this low-level header must only be included from ylapack.h
#endif

_YLPK_BEGIN_DECLS

/*---------------------------------------------------------------------------*/
/* BLAS LEVEL 0 */

#if (FORTRAN_STYLE != 1)
# define XERBLA   BLAS_NAME(xerbla,XERBLA)
# define LSAME    BLAS_NAME(lsame,LSAME)
#endif /* FORTRAN_STYLE */

extern VOID XERBLA(CONST CHARACTER *srcname, CONST INTEGER *info,
                   CONST int len_srcname);
extern LOGICAL LSAME(CONST CHARACTER *ca, CONST CHARACTER *cb,
                     CONST int len_ca, CONST int len_cb);

/*---------------------------------------------------------------------------*/
/* BLAS LEVEL 1 */

/* _ROTG  -- Generate plane rotation.
 * _ROTMG -- Generate modified plane rotation.
 * _ROT   -- Apply plane rotation.
 * _ROTM  -- Apply modified plane rotation.
 * _SWAP  -- Exchange contents of 2 vectors.
 * _SCAL  -- Scale a vector by a scalar.
 * _COPY  -- Copy contents of a vector.
 * _AXPY  -- Y := ALPHA*X + Y
 * _DOT  -- Compute dot product of 2 vectors: x^T*y
 * _DOTU -- Compute dot product of 2 vectors: x^T*y
 * _DOTC -- Compute dot product of 2 vectors: x^H*y
 * __DOT -- alpha + x^T*y
 * _NRM2 -- Euclidean (L-2) norm of a vector.
 * _ASUM -- L-1 norm of a vector: ||Re(x)||_1 + ||Im(x)||_1
 * I_AMAX -- Index of 1st element with maximum |Re(x_k)| + |Im(x_k)|
 */

#if (FORTRAN_STYLE != 1)
# define SCABS1   BLAS_NAME(scabs1,SCABS1)
# define SROTG    BLAS_NAME(srotg,SROTG)
# define SROTM    BLAS_NAME(srotm,SROTM)
# define SROTMG   BLAS_NAME(srotmg,SROTMG)
# define SROT     BLAS_NAME(srot,SROT)
# define SCOPY    BLAS_NAME(scopy,SCOPY)
# define SSWAP    BLAS_NAME(sswap,SSWAP)
# define SSCAL    BLAS_NAME(sscal,SSCAL)
# define SAXPY    BLAS_NAME(saxpy,SAXPY)
# define SDOT     BLAS_NAME(sdot,SDOT)
# define SDSDOT   BLAS_NAME(sdsdot,SDSDOT)
# define DSDOT    BLAS_NAME(dsdot,DSDOT)
# define SNRM2    BLAS_NAME(snrm2,SNRM2)
# define SASUM    BLAS_NAME(sasum,SASUM)
# define ISAMAX   BLAS_NAME(isamax,ISAMAX)
#endif /* FORTRAN_STYLE */

extern float   SCABS1(CONST COMPLEX8_PTR c);
extern VOID     SROTG(float *a, float *b, float *c, float *s);
extern VOID     SROTM(CONST INTEGER *n, float *x, CONST INTEGER *incx,
                      float *y, CONST INTEGER *incy, CONST float param[5]);
extern VOID    SROTMG(float *d1, float *d2, float *x1, CONST float *y1,
                      float param[5]);
extern VOID      SROT(CONST INTEGER *n, float *x, CONST INTEGER *incx, float *y,
                      CONST INTEGER *incy, CONST float *c, CONST float *s);
extern VOID     SCOPY(CONST INTEGER *n, CONST float *x, CONST INTEGER *incx,
                      float *y, CONST INTEGER *incy);
extern VOID     SSWAP(CONST INTEGER *n, float *x, CONST INTEGER *incx, float *y,
                      CONST INTEGER *incy);
extern VOID     SSCAL(CONST INTEGER *n, CONST float *a, float *x,
                      CONST INTEGER *incx);
extern VOID     SAXPY(CONST INTEGER *n, CONST float *alpha, CONST float *x,
                      CONST INTEGER *incx, float *y, CONST INTEGER *incy);
extern float     SDOT(CONST INTEGER *n, CONST float *x, CONST INTEGER *incx,
                      CONST float *y, CONST INTEGER *incy);
extern float   SDSDOT(CONST INTEGER *n, CONST float *sb, CONST float *x,
                      CONST INTEGER *incx, CONST float *y,
                      CONST INTEGER *incy);
extern double   DSDOT(CONST INTEGER *n, CONST float *x, CONST INTEGER *incx,
                      CONST float *y, CONST INTEGER *incy);
extern float    SNRM2(CONST INTEGER *n, CONST float *x, CONST INTEGER *incx);
extern float    SASUM(CONST INTEGER *n, CONST float *x, CONST INTEGER *incx);
extern INTEGER ISAMAX(CONST INTEGER *n, CONST float *x, CONST INTEGER *incx);

#if (FORTRAN_STYLE != 1)
# define DCABS1   BLAS_NAME(dcabs1,DCABS1)
# define DROTG    BLAS_NAME(drotg,DROTG)
# define DROTM    BLAS_NAME(drotm,DROTM)
# define DROTMG   BLAS_NAME(drotmg,DROTMG)
# define DROT     BLAS_NAME(drot,DROT)
# define DCOPY    BLAS_NAME(dcopy,DCOPY)
# define DSWAP    BLAS_NAME(dswap,DSWAP)
# define DSCAL    BLAS_NAME(dscal,DSCAL)
# define DAXPY    BLAS_NAME(daxpy,DAXPY)
# define DDOT     BLAS_NAME(ddot,DDOT)
# define DNRM2    BLAS_NAME(dnrm2,DNRM2)
# define DASUM    BLAS_NAME(dasum,DASUM)
# define IDAMAX   BLAS_NAME(idamax,IDAMAX)
#endif /* FORTRAN_STYLE */

extern double   DCABS1(CONST COMPLEX16_PTR z);
extern VOID      DROTG(double *a, double *b, double *c, double *s);
extern VOID      DROTM(CONST INTEGER *n, double *x, CONST INTEGER *incx,
                       double *y, CONST INTEGER *incy, CONST double param[5]);
extern VOID     DROTMG(double *d1, double *d2, double *x1, CONST double *y1,
                       double *param[5]);
extern VOID       DROT(CONST INTEGER *n, double *x, CONST INTEGER *incx,
                       double *y, CONST INTEGER *incy, CONST double *c,
                       CONST double *s);
extern VOID      DCOPY(CONST INTEGER *n, CONST double *x, CONST INTEGER *incx,
                       double *y, CONST INTEGER *incy);
extern VOID      DSWAP(CONST INTEGER *n, double *x, CONST INTEGER *incx,
                       double *y, CONST INTEGER *incy);
extern VOID      DSCAL(CONST INTEGER *n, CONST double *a, double *x,
                       CONST INTEGER *incx);
extern VOID      DAXPY(CONST INTEGER *n, CONST double *alpha, CONST double *x,
                       CONST INTEGER *incx, double *y, CONST INTEGER *incy);
extern double     DDOT(CONST INTEGER *n, CONST double *x, CONST INTEGER *incx,
                       CONST double *y, CONST INTEGER *incy);
extern double    DNRM2(CONST INTEGER *n, CONST double *x, CONST INTEGER *incx);
extern double    DASUM(CONST INTEGER *n, CONST double *x, CONST INTEGER *incx);
extern INTEGER  IDAMAX(CONST INTEGER *n, CONST double *x, CONST INTEGER *incx);

#define SCALAR(TYPE, NAME) TYPE##_PTR NAME
#define VECTOR(TYPE, NAME) TYPE##_PTR NAME, CONST INTEGER *inc##NAME
#define MATRIX(TYPE, NAME) TYPE##_PTR NAME, CONST INTEGER *ld##NAME
#define float_PTR     float  *
#define double_PTR    double *
#if 0
# define complex8_PTR  float  *
# define complex16_PTR double *
#else
# define complex8_PTR  void *
# define complex16_PTR void *
#endif
extern float     SNRM2(CONST INTEGER *n, CONST VECTOR(float,  x));
extern double    DNRM2(CONST INTEGER *n, CONST VECTOR(double, x));
extern float     SASUM(CONST INTEGER *n, CONST VECTOR(float,  x));
extern double    DASUM(CONST INTEGER *n, CONST VECTOR(double, x));

#if (FORTRAN_STYLE != 1)
# define CROTG    BLAS_NAME(crotg,CROTG)
# define CSROT    BLAS_NAME(csrot,CSROT)
# define CROT     BLAS_NAME(crot,CROT)
# define CCOPY    BLAS_NAME(ccopy,CCOPY)
# define CSWAP    BLAS_NAME(cswap,CSWAP)
# define CSSCAL   BLAS_NAME(csscal,CSSCAL)
# define CSCAL    BLAS_NAME(cscal,CSCAL)
# define CAXPY    BLAS_NAME(caxpy,CAXPY)
# define CDOTC    BLAS_NAME(cdotc,CDOTC)
# define CDOTU    BLAS_NAME(cdotu,CDOTU)
# define SCNRM2   BLAS_NAME(scnrm2,SCNRM2)
# define SCASUM   BLAS_NAME(scasum,SCASUM)
# define ICAMAX   BLAS_NAME(icamax,ICAMAX)
#endif /* FORTRAN_STYLE */

extern VOID     CROTG(COMPLEX8_PTR a, CONST COMPLEX8_PTR b, float *c,
                      COMPLEX8_PTR s);
extern VOID     CSROT(CONST INTEGER *n, COMPLEX8_PTR x, CONST INTEGER *incx,
                      COMPLEX8_PTR y, CONST INTEGER *incy, CONST float *c,
                      CONST float *s);
extern VOID      CROT(CONST INTEGER *n, COMPLEX8_PTR x, CONST INTEGER *incx,
                      COMPLEX8_PTR y, CONST INTEGER *incy, CONST float *c,
                      CONST COMPLEX8_PTR s);
extern VOID     CCOPY(CONST INTEGER *n, CONST COMPLEX8_PTR x, CONST INTEGER *incx,
                      COMPLEX8_PTR y, CONST INTEGER *incy);
extern VOID     CSWAP(CONST INTEGER *n, COMPLEX8_PTR x, CONST INTEGER *incx,
                      COMPLEX8_PTR y, CONST INTEGER *incy);
extern VOID    CSSCAL(CONST INTEGER *n, CONST float *a, COMPLEX8_PTR x,
                      CONST INTEGER *incx);
extern VOID     CSCAL(CONST INTEGER *n, CONST COMPLEX8_PTR a, COMPLEX8_PTR x,
                      CONST INTEGER *incx);
extern VOID     CAXPY(CONST INTEGER *n, CONST COMPLEX8_PTR alpha,
                      CONST COMPLEX8_PTR x, CONST INTEGER *incx, COMPLEX8_PTR y,
                      CONST INTEGER *incy);
extern VOID     CDOTC(COMPLEX8_PTR result, CONST INTEGER *n,
                      CONST COMPLEX8_PTR x, CONST INTEGER *incx,
                      CONST COMPLEX8_PTR y, CONST INTEGER *incy);
extern VOID     CDOTU(COMPLEX8_PTR result, CONST INTEGER *n,
                      CONST COMPLEX8_PTR x, CONST INTEGER *incx,
                      CONST COMPLEX8_PTR y, CONST INTEGER *incy);
extern float   SCNRM2(CONST INTEGER *n, CONST COMPLEX8_PTR x, CONST INTEGER *incx);
extern float   SCASUM(CONST INTEGER *n, CONST COMPLEX8_PTR x, CONST INTEGER *incx);
extern INTEGER ICAMAX(CONST INTEGER *n, CONST COMPLEX8_PTR x,
                      CONST INTEGER *incx);

#if (FORTRAN_STYLE != 1)
# define ZROTG    BLAS_NAME(zrotg,ZROTG)
# define ZDROT    BLAS_NAME(zdrot,ZDROT)
# define ZROT     BLAS_NAME(zrot,ZROT)
# define ZCOPY    BLAS_NAME(zcopy,ZCOPY)
# define ZSWAP    BLAS_NAME(zswap,ZSWAP)
# define ZDSCAL   BLAS_NAME(zdscal,ZDSCAL)
# define ZSCAL    BLAS_NAME(zscal,ZSCAL)
# define ZAXPY    BLAS_NAME(zaxpy,ZAXPY)
# define ZDOTC    BLAS_NAME(zdotc,ZDOTC)
# define ZDOTU    BLAS_NAME(zdotu,ZDOTU)
# define DZNRM2   BLAS_NAME(dznrm2,DZNRM2)
# define DZASUM   BLAS_NAME(dzasum,DZASUM)
# define IZAMAX   BLAS_NAME(izamax,IZAMAX)
#endif /* FORTRAN_STYLE */

extern VOID     ZROTG(COMPLEX16_PTR a, CONST COMPLEX16_PTR b, double *c,
                      COMPLEX16_PTR s);
extern VOID     ZDROT(CONST INTEGER *n, COMPLEX16_PTR x, CONST INTEGER *incx,
                      COMPLEX16_PTR y, CONST INTEGER *incy, CONST double *c,
                      CONST double *s);
extern VOID      ZROT(CONST INTEGER *n, COMPLEX16_PTR x, CONST INTEGER *incx,
                      COMPLEX16_PTR y, CONST INTEGER *incy, CONST double *c,
                      CONST COMPLEX16_PTR s);
extern VOID     ZCOPY(CONST INTEGER *n,
                      CONST COMPLEX16_PTR x, CONST INTEGER *incx,
                      COMPLEX16_PTR y, CONST INTEGER *incy);
extern VOID     ZSWAP(CONST INTEGER *n, COMPLEX16_PTR x, CONST INTEGER *incx,
                      COMPLEX16_PTR y, CONST INTEGER *incy);
extern VOID    ZDSCAL(CONST INTEGER *n, CONST double *a, COMPLEX16_PTR x,
                      CONST INTEGER *incx);
extern VOID     ZSCAL(CONST INTEGER *n, CONST COMPLEX16_PTR a, COMPLEX16_PTR x,
                      CONST INTEGER *incx);
extern VOID     ZAXPY(CONST INTEGER *n, CONST COMPLEX16_PTR alpha,
                      CONST COMPLEX16_PTR x, CONST INTEGER *incx, COMPLEX16_PTR y,
                      CONST INTEGER *incy);
extern VOID     ZDOTC(COMPLEX16_PTR result, CONST INTEGER *n,
                      CONST COMPLEX16_PTR x, CONST INTEGER *incx,
                      CONST COMPLEX16_PTR y, CONST INTEGER *incy);
extern VOID     ZDOTU(COMPLEX16_PTR result, CONST INTEGER *n,
                      CONST COMPLEX16_PTR x, CONST INTEGER *incx,
                      CONST COMPLEX16_PTR y, CONST INTEGER *incy);
extern double  DZNRM2(CONST INTEGER *n, CONST COMPLEX16_PTR x,
                      CONST INTEGER *incx);
extern double  DZASUM(CONST INTEGER *n, CONST COMPLEX16_PTR x,
                      CONST INTEGER *incx);
extern INTEGER IZAMAX(CONST INTEGER *n, CONST COMPLEX16_PTR x,
                      CONST INTEGER *incx);

/*---------------------------------------------------------------------------*/
/* BLAS LEVEL 2 */

#if (FORTRAN_STYLE != 1)
# define SGEMV   BLAS_NAME(sgemv,SGEMV)
# define SGBMV   BLAS_NAME(sgbmv,SGBMV)
# define SSYMV   BLAS_NAME(ssymv,SSYMV)
# define SSBMV   BLAS_NAME(ssbmv,SSBMV)
# define SSPMV   BLAS_NAME(sspmv,SSPMV)
# define STRMV   BLAS_NAME(strmv,STRMV)
# define STBMV   BLAS_NAME(stbmv,STBMV)
# define STPMV   BLAS_NAME(stpmv,STPMV)
# define STRSV   BLAS_NAME(strsv,STRSV)
# define STBSV   BLAS_NAME(stbsv,STBSV)
# define STPSV   BLAS_NAME(stpsv,STPSV)
# define SGER    BLAS_NAME(sger,SGER)
# define SSYR    BLAS_NAME(ssyr,SSYR)
# define SSPR    BLAS_NAME(sspr,SSPR)
# define SSYR2   BLAS_NAME(ssyr2,SSYR2)
# define SSPR2   BLAS_NAME(sspr2,SSPR2)
#endif /* FORTRAN_STYLE */

extern VOID SGBMV(CONST CHARACTER *trans, CONST INTEGER *m,
                  CONST INTEGER *n, CONST INTEGER *kl, CONST INTEGER *ku,
                  CONST float *alpha, CONST float *a, CONST INTEGER *lda,
                  CONST float *x, CONST INTEGER *incx, CONST float *beta,
                  float *y, CONST INTEGER *incy);
extern VOID SGEMV(CONST CHARACTER *trans, CONST INTEGER *m,
                  CONST INTEGER *n, CONST float *alpha, CONST float *a,
                  CONST INTEGER *lda, CONST float *x, CONST INTEGER *incx,
                  CONST float *beta, float *y, CONST INTEGER *incy);
extern VOID  SGER(CONST INTEGER *m, CONST INTEGER *n, CONST float *alpha,
                  CONST float *x, CONST INTEGER *incx, CONST float *y,
                  CONST INTEGER *incy, float *a, CONST INTEGER *lda);
extern VOID SSBMV(CONST CHARACTER *uplo, CONST INTEGER *n, CONST INTEGER *k,
                  CONST float *alpha, CONST float *a, CONST INTEGER *lda,
                  CONST float *x, CONST INTEGER *incx, CONST float *beta,
                  float *y, CONST INTEGER *incy);
extern VOID SSPMV(CONST CHARACTER *uplo, CONST INTEGER *n,
                  CONST float *alpha, CONST float *ap, CONST float *x,
                  CONST INTEGER *incx, CONST float *beta, float *y,
                  CONST INTEGER *incy);
extern VOID  SSPR(CONST CHARACTER *uplo, CONST INTEGER *n, CONST float *alpha,
                  CONST float *x, CONST INTEGER *incx, float *ap);
extern VOID SSPR2(CONST CHARACTER *uplo, CONST INTEGER *n,
                  CONST float *alpha, CONST float *x, CONST INTEGER *incx,
                  CONST float *y, CONST INTEGER *incy, float *ap);
extern VOID SSYMV(CONST CHARACTER *uplo, CONST INTEGER *n,
                  CONST float *alpha, CONST float *a, CONST INTEGER *lda,
                  CONST float *x, CONST INTEGER *incx, CONST float *beta,
                  float *y, CONST INTEGER *incy);
extern VOID  SSYR(CONST CHARACTER *uplo, CONST INTEGER *n, CONST float *alpha,
                  CONST float *x, CONST INTEGER *incx, float *a,
                  CONST INTEGER *lda);
extern VOID SSYR2(CONST CHARACTER *uplo, CONST INTEGER *n,
                  CONST float *alpha, CONST float *x, CONST INTEGER *incx,
                  CONST float *y, CONST INTEGER *incy, float *a,
                  CONST INTEGER *lda);
extern VOID STBMV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                  CONST CHARACTER *diag, CONST INTEGER *n, CONST INTEGER *k,
                  CONST float *a, CONST INTEGER *lda, float *x,
                  CONST INTEGER *incx);
extern VOID STBSV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                  CONST CHARACTER *diag, CONST INTEGER *n, CONST INTEGER *k,
                  CONST float *a, CONST INTEGER *lda, float *x,
                  CONST INTEGER *incx);
extern VOID STPMV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                  CONST CHARACTER *diag, CONST INTEGER *n, CONST float *ap,
                  float *x, CONST INTEGER *incx);
extern VOID STPSV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                  CONST CHARACTER *diag, CONST INTEGER *n, CONST float *ap,
                  float *x, CONST INTEGER *incx);
extern VOID STRMV(CONST CHARACTER *uplo, CONST CHARACTER *transa,
                  CONST CHARACTER *diag, CONST INTEGER *n, CONST float *a,
                  CONST INTEGER *lda, float *b, CONST INTEGER *incx);
extern VOID STRSV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                  CONST CHARACTER *diag, CONST INTEGER *n, CONST float *a,
                  CONST INTEGER *lda, float *x, CONST INTEGER *incx);

#if (FORTRAN_STYLE != 1)
# define DGEMV   BLAS_NAME(dgemv,DGEMV)
# define DGBMV   BLAS_NAME(dgbmv,DGBMV)
# define DSYMV   BLAS_NAME(dsymv,DSYMV)
# define DSBMV   BLAS_NAME(dsbmv,DSBMV)
# define DSPMV   BLAS_NAME(dspmv,DSPMV)
# define DTRMV   BLAS_NAME(dtrmv,DTRMV)
# define DTBMV   BLAS_NAME(dtbmv,DTBMV)
# define DTPMV   BLAS_NAME(dtpmv,DTPMV)
# define DTRSV   BLAS_NAME(dtrsv,DTRSV)
# define DTBSV   BLAS_NAME(dtbsv,DTBSV)
# define DTPSV   BLAS_NAME(dtpsv,DTPSV)
# define DGER    BLAS_NAME(dger,DGER)
# define DSYR    BLAS_NAME(dsyr,DSYR)
# define DSPR    BLAS_NAME(dspr,DSPR)
# define DSYR2   BLAS_NAME(dsyr2,DSYR2)
# define DSPR2   BLAS_NAME(dspr2,DSPR2)
#endif /* FORTRAN_STYLE */

extern VOID DGBMV(CONST CHARACTER *trans, CONST INTEGER *m,
                  CONST INTEGER *n, CONST INTEGER *kl, CONST INTEGER *ku,
                  CONST double *alpha, CONST double *a, CONST INTEGER *lda,
                  CONST double *x, CONST INTEGER *incx, CONST double *beta,
                  double *y, CONST INTEGER *incy);
extern VOID DGEMV(CONST CHARACTER *trans, CONST INTEGER *m,
                  CONST INTEGER *n, CONST double *alpha, CONST double *a,
                  CONST INTEGER *lda, CONST double *x, CONST INTEGER *incx,
                  CONST double *beta, double *y, CONST INTEGER *incy);
extern VOID  DGER(CONST INTEGER *m, CONST INTEGER *n, CONST double *alpha,
                  CONST double *x, CONST INTEGER *incx, CONST double *y,
                  CONST INTEGER *incy, double *a, CONST INTEGER *lda);
extern VOID DSBMV(CONST CHARACTER *uplo, CONST INTEGER *n, CONST INTEGER *k,
                  CONST double *alpha, CONST double *a, CONST INTEGER *lda,
                  CONST double *x, CONST INTEGER *incx, CONST double *beta,
                  double *y, CONST INTEGER *incy);
extern VOID DSPMV(CONST CHARACTER *uplo, CONST INTEGER *n,
                  CONST double *alpha, CONST double *ap, CONST double *x,
                  CONST INTEGER *incx, CONST double *beta, double *y,
                  CONST INTEGER *incy);
extern VOID  DSPR(CONST CHARACTER *uplo, CONST INTEGER *n,
                  CONST double *alpha, CONST double *x, CONST INTEGER *incx,
                  double *ap);
extern VOID DSPR2(CONST CHARACTER *uplo, CONST INTEGER *n,
                  CONST double *alpha, CONST double *x, CONST INTEGER *incx,
                  CONST double *y, CONST INTEGER *incy, double *ap);
extern VOID DSYMV(CONST CHARACTER *uplo, CONST INTEGER *n,
                  CONST double *alpha, CONST double *a, CONST INTEGER *lda,
                  CONST double *x, CONST INTEGER *incx, CONST double *beta,
                  double *y, CONST INTEGER *incy);
extern VOID  DSYR(CONST CHARACTER *uplo, CONST INTEGER *n,
                  CONST double *alpha, CONST double *x, CONST INTEGER *incx,
                  double *a, CONST INTEGER *lda);
extern VOID DSYR2(CONST CHARACTER *uplo, CONST INTEGER *n,
                  CONST double *alpha, CONST double *x, CONST INTEGER *incx,
                  CONST double *y, CONST INTEGER *incy, double *a,
                  CONST INTEGER *lda);
extern VOID DTBMV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                  CONST CHARACTER *diag, CONST INTEGER *n, CONST INTEGER *k,
                  CONST double *a, CONST INTEGER *lda, double *x,
                  CONST INTEGER *incx);
extern VOID DTBSV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                  CONST CHARACTER *diag, CONST INTEGER *n, CONST INTEGER *k,
                  CONST double *a, CONST INTEGER *lda, double *x,
                  CONST INTEGER *incx);
extern VOID DTPMV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                  CONST CHARACTER *diag, CONST INTEGER *n, CONST double *ap,
                  double *x, CONST INTEGER *incx);
extern VOID DTPSV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                  CONST CHARACTER *diag, CONST INTEGER *n, CONST double *ap,
                  double *x, CONST INTEGER *incx);
extern VOID DTRMV(CONST CHARACTER *uplo, CONST CHARACTER *transa,
                  CONST CHARACTER *diag, CONST INTEGER *n, CONST double *a,
                  CONST INTEGER *lda, double *b, CONST INTEGER *incx);
extern VOID DTRSV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                  CONST CHARACTER *diag, CONST INTEGER *n, CONST double *a,
                  CONST INTEGER *lda, double *x, CONST INTEGER *incx);

#if (FORTRAN_STYLE != 1)
# define CGBMV   BLAS_NAME(cgbmv,CGBMV)
# define CGEMV   BLAS_NAME(cgemv,CGEMV)
# define CGERC   BLAS_NAME(cgerc,CGERC)
# define CGERU   BLAS_NAME(cgeru,CGERU)
# define CHBMV   BLAS_NAME(chbmv,CHBMV)
# define CHEMV   BLAS_NAME(chemv,CHEMV)
# define CHER    BLAS_NAME(cher,CHER)
# define CHER2   BLAS_NAME(cher2,CHER2)
# define CHPMV   BLAS_NAME(chpmv,CHPMV)
# define CHPR    BLAS_NAME(chpr,CHPR)
# define CHPR2   BLAS_NAME(chpr2,CHPR2)
# define CTBMV   BLAS_NAME(ctbmv,CTBMV)
# define CTBSV   BLAS_NAME(ctbsv,CTBSV)
# define CTPMV   BLAS_NAME(ctpmv,CTPMV)
# define CTPSV   BLAS_NAME(ctpsv,CTPSV)
# define CTRMV   BLAS_NAME(ctrmv,CTRMV)
# define CTRSV   BLAS_NAME(ctrsv,CTRSV)
# define SCGEMV  BLAS_NAME(scgemv,SCGEMV)
#endif /* FORTRAN_STYLE */

extern VOID   CGBMV(CONST CHARACTER *trans, CONST INTEGER *m,
                    CONST INTEGER *n, CONST INTEGER *kl, CONST INTEGER *ku,
                    CONST COMPLEX8_PTR alpha, CONST COMPLEX8_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX8_PTR x,
                    CONST INTEGER *incx, CONST COMPLEX8_PTR beta, COMPLEX8_PTR y,
                    CONST INTEGER *incy);
extern VOID   CGEMV(CONST CHARACTER *trans, CONST INTEGER *m,
                    CONST INTEGER *n, CONST COMPLEX8_PTR alpha,
                    CONST COMPLEX8_PTR a, CONST INTEGER *lda, CONST COMPLEX8_PTR x,
                    CONST INTEGER *incx, CONST COMPLEX8_PTR beta, COMPLEX8_PTR y,
                    CONST INTEGER *incy);
extern VOID   CGERC(CONST INTEGER *m, CONST INTEGER *n, CONST COMPLEX8_PTR alpha,
                    CONST COMPLEX8_PTR x, CONST INTEGER *incx, CONST COMPLEX8_PTR y,
                    CONST INTEGER *incy, COMPLEX8_PTR a, CONST INTEGER *lda);
extern VOID   CGERU(CONST INTEGER *m, CONST INTEGER *n, CONST COMPLEX8_PTR alpha,
                    CONST COMPLEX8_PTR x, CONST INTEGER *incx, CONST COMPLEX8_PTR y,
                    CONST INTEGER *incy, COMPLEX8_PTR a, CONST INTEGER *lda);
extern VOID   CHBMV(CONST CHARACTER *uplo, CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX8_PTR alpha, CONST COMPLEX8_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX8_PTR x,
                    CONST INTEGER *incx, CONST COMPLEX8_PTR beta, COMPLEX8_PTR y,
                    CONST INTEGER *incy);
extern VOID   CHEMV(CONST CHARACTER *uplo, CONST INTEGER *n,
                    CONST COMPLEX8_PTR alpha, CONST COMPLEX8_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX8_PTR x,
                    CONST INTEGER *incx, CONST COMPLEX8_PTR beta, COMPLEX8_PTR y,
                    CONST INTEGER *incy);
extern VOID    CHER(CONST CHARACTER *uplo, CONST INTEGER *n, CONST float *alpha,
                    CONST COMPLEX8_PTR x, CONST INTEGER *incx, COMPLEX8_PTR a,
                    CONST INTEGER *lda);
extern VOID   CHER2(CONST CHARACTER *uplo, CONST INTEGER *n,
                    CONST COMPLEX8_PTR alpha, CONST COMPLEX8_PTR x,
                    CONST INTEGER *incx, CONST COMPLEX8_PTR y,
                    CONST INTEGER *incy, COMPLEX8_PTR a, CONST INTEGER *lda);
extern VOID   CHPMV(CONST CHARACTER *uplo, CONST INTEGER *n,
                    CONST COMPLEX8_PTR alpha, CONST COMPLEX8_PTR ap,
                    CONST COMPLEX8_PTR x, CONST INTEGER *incx,
                    CONST COMPLEX8_PTR beta, COMPLEX8_PTR y, CONST INTEGER *incy);
extern VOID    CHPR(CONST CHARACTER *uplo, CONST INTEGER *n, CONST float *alpha,
                    CONST COMPLEX8_PTR x, CONST INTEGER *incx, COMPLEX8_PTR ap);
extern VOID   CHPR2(CONST CHARACTER *uplo, CONST INTEGER *n,
                    CONST COMPLEX8_PTR alpha, CONST COMPLEX8_PTR x,
                    CONST INTEGER *incx, CONST COMPLEX8_PTR y,
                    CONST INTEGER *incy, COMPLEX8_PTR ap);
extern VOID   CTBMV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST CHARACTER *diag, CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX8_PTR a, CONST INTEGER *lda, COMPLEX8_PTR x,
                    CONST INTEGER *incx);
extern VOID   CTBSV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST CHARACTER *diag, CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX8_PTR a, CONST INTEGER *lda, COMPLEX8_PTR x,
                    CONST INTEGER *incx);
extern VOID   CTPMV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST CHARACTER *diag, CONST INTEGER *n,
                    CONST COMPLEX8_PTR ap, COMPLEX8_PTR x, CONST INTEGER *incx);
extern VOID   CTPSV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST CHARACTER *diag, CONST INTEGER *n,
                    CONST COMPLEX8_PTR ap, COMPLEX8_PTR x, CONST INTEGER *incx);
extern VOID   CTRMV(CONST CHARACTER *uplo, CONST CHARACTER *transa,
                    CONST CHARACTER *diag, CONST INTEGER *n,
                    CONST COMPLEX8_PTR a, CONST INTEGER *lda, COMPLEX8_PTR b,
                    CONST INTEGER *incx);
extern VOID   CTRSV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST CHARACTER *diag, CONST INTEGER *n,
                    CONST COMPLEX8_PTR a, CONST INTEGER *lda, COMPLEX8_PTR x,
                    CONST INTEGER *incx);
extern VOID  SCGEMV(CONST CHARACTER *trans, CONST INTEGER *m,
                    CONST INTEGER *n, CONST COMPLEX8_PTR alpha, CONST float *a,
                    CONST INTEGER *lda, CONST COMPLEX8_PTR x,
                    CONST INTEGER *incx, CONST COMPLEX8_PTR beta, COMPLEX8_PTR y,
                    CONST INTEGER *incy);

#if (FORTRAN_STYLE != 1)
# define ZGBMV   BLAS_NAME(zgbmv,ZGBMV)
# define ZGEMV   BLAS_NAME(zgemv,ZGEMV)
# define ZGERC   BLAS_NAME(zgerc,ZGERC)
# define ZGERU   BLAS_NAME(zgeru,ZGERU)
# define ZHBMV   BLAS_NAME(zhbmv,ZHBMV)
# define ZHEMV   BLAS_NAME(zhemv,ZHEMV)
# define ZHER    BLAS_NAME(zher,ZHER)
# define ZHER2   BLAS_NAME(zher2,ZHER2)
# define ZHPMV   BLAS_NAME(zhpmv,ZHPMV)
# define ZHPR    BLAS_NAME(zhpr,ZHPR)
# define ZHPR2   BLAS_NAME(zhpr2,ZHPR2)
# define ZTBMV   BLAS_NAME(ztbmv,ZTBMV)
# define ZTBSV   BLAS_NAME(ztbsv,ZTBSV)
# define ZTPMV   BLAS_NAME(ztpmv,ZTPMV)
# define ZTPSV   BLAS_NAME(ztpsv,ZTPSV)
# define ZTRMV   BLAS_NAME(ztrmv,ZTRMV)
# define ZTRSV   BLAS_NAME(ztrsv,ZTRSV)
# define DZGEMV  BLAS_NAME(dzgemv,DZGEMV)
#endif /* FORTRAN_STYLE */

extern VOID   ZGBMV(CONST CHARACTER *trans, CONST INTEGER *m,
                    CONST INTEGER *n, CONST INTEGER *kl, CONST INTEGER *ku,
                    CONST COMPLEX16_PTR alpha,
                    CONST COMPLEX16_PTR a, CONST INTEGER *lda,
                    CONST COMPLEX16_PTR x, CONST INTEGER *incx,
                    CONST COMPLEX16_PTR beta,
                    COMPLEX16_PTR y, CONST INTEGER *incy);
extern VOID   ZGEMV(CONST CHARACTER *trans, CONST INTEGER *m,
                    CONST INTEGER *n, CONST COMPLEX16_PTR alpha,
                    CONST COMPLEX16_PTR a, CONST INTEGER *lda,
                    CONST COMPLEX16_PTR x, CONST INTEGER *incx,
                    CONST COMPLEX16_PTR beta, COMPLEX16_PTR y, CONST INTEGER *incy);
extern VOID   ZGERC(CONST INTEGER *m, CONST INTEGER *n,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR x,
                    CONST INTEGER *incx, CONST COMPLEX16_PTR y,
                    CONST INTEGER *incy, COMPLEX16_PTR a, CONST INTEGER *lda);
extern VOID   ZGERU(CONST INTEGER *m, CONST INTEGER *n,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR x,
                    CONST INTEGER *incx, CONST COMPLEX16_PTR y,
                    CONST INTEGER *incy, COMPLEX16_PTR a, CONST INTEGER *lda);
extern VOID   ZHBMV(CONST CHARACTER *uplo, CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX16_PTR x,
                    CONST INTEGER *incx, CONST COMPLEX16_PTR beta, COMPLEX16_PTR y,
                    CONST INTEGER *incy);
extern VOID   ZHEMV(CONST CHARACTER *uplo, CONST INTEGER *n,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX16_PTR x,
                    CONST INTEGER *incx, CONST COMPLEX16_PTR beta, COMPLEX16_PTR y,
                    CONST INTEGER *incy);
extern VOID    ZHER(CONST CHARACTER *uplo, CONST INTEGER *n,
                    CONST double *alpha, CONST COMPLEX16_PTR x,
                    CONST INTEGER *incx, COMPLEX16_PTR a, CONST INTEGER *lda);
extern VOID   ZHER2(CONST CHARACTER *uplo, CONST INTEGER *n,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR x,
                    CONST INTEGER *incx, CONST COMPLEX16_PTR y,
                    CONST INTEGER *incy, COMPLEX16_PTR a, CONST INTEGER *lda);
extern VOID   ZHPMV(CONST CHARACTER *uplo, CONST INTEGER *n,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR ap,
                    CONST COMPLEX16_PTR x, CONST INTEGER *incx,
                    CONST COMPLEX16_PTR beta, COMPLEX16_PTR y, CONST INTEGER *incy);
extern VOID    ZHPR(CONST CHARACTER *uplo, CONST INTEGER *n,
                    CONST double *alpha, CONST COMPLEX16_PTR x,
                    CONST INTEGER *incx, COMPLEX16_PTR ap);
extern VOID   ZHPR2(CONST CHARACTER *uplo, CONST INTEGER *n,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR x,
                    CONST INTEGER *incx, CONST COMPLEX16_PTR y,
                    CONST INTEGER *incy, COMPLEX16_PTR ap);
extern VOID   ZTBMV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST CHARACTER *diag, CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX16_PTR a, CONST INTEGER *lda, COMPLEX16_PTR x,
                    CONST INTEGER *incx);
extern VOID   ZTBSV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST CHARACTER *diag, CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX16_PTR a, CONST INTEGER *lda, COMPLEX16_PTR x,
                    CONST INTEGER *incx);
extern VOID   ZTPMV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST CHARACTER *diag, CONST INTEGER *n,
                    CONST COMPLEX16_PTR ap, COMPLEX16_PTR x,
                    CONST INTEGER *incx);
extern VOID   ZTPSV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST CHARACTER *diag, CONST INTEGER *n, COMPLEX16_PTR ap,
                    COMPLEX16_PTR x, CONST INTEGER *incx);
extern VOID   ZTRMV(CONST CHARACTER *uplo, CONST CHARACTER *transa,
                    CONST CHARACTER *diag, CONST INTEGER *n,
                    CONST COMPLEX16_PTR a, CONST INTEGER *lda, COMPLEX16_PTR b,
                    CONST INTEGER *incx);
extern VOID   ZTRSV(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST CHARACTER *diag, CONST INTEGER *n,
                    CONST COMPLEX16_PTR a, CONST INTEGER *lda, COMPLEX16_PTR x,
                    CONST INTEGER *incx);
extern VOID  DZGEMV(CONST CHARACTER *trans, CONST INTEGER *m,
                    CONST INTEGER *n, CONST COMPLEX16_PTR alpha,
                    CONST double *a, CONST INTEGER *lda,
                    CONST COMPLEX16_PTR x, CONST INTEGER *incx,
                    CONST COMPLEX16_PTR beta, COMPLEX16_PTR y,
                    CONST INTEGER *incy);

/*---------------------------------------------------------------------------*/
/* BLAS LEVEL 3 */

#if (FORTRAN_STYLE != 1)
# define SGEMM   BLAS_NAME(sgemm,SGEMM)
# define SSYMM   BLAS_NAME(ssymm,SSYMM)
# define SSYRK   BLAS_NAME(ssyrk,SSYRK)
# define SSYR2K  BLAS_NAME(ssyr2k,SSYR2K)
# define STRMM   BLAS_NAME(strmm,STRMM)
# define STRSM   BLAS_NAME(strsm,STRSM)
#endif /* FORTRAN_STYLE */

extern VOID  SGEMM(CONST CHARACTER *transa, CONST CHARACTER *transb,
                   CONST INTEGER *m, CONST INTEGER *n, CONST INTEGER *k,
                   CONST float *alpha, CONST float *a, CONST INTEGER *lda,
                   CONST float *b, CONST INTEGER *ldb, CONST float *beta,
                   float *c, CONST INTEGER *ldc);
extern VOID  SSYMM(CONST CHARACTER *side, CONST CHARACTER *uplo,
                   CONST INTEGER *m, CONST INTEGER *n, CONST float *alpha,
                   CONST float *a, CONST INTEGER *lda, CONST float *b,
                   CONST INTEGER *ldb, CONST float *beta, float *c,
                   CONST INTEGER *ldc);
extern VOID  SSYRK(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                   CONST INTEGER *n, CONST INTEGER *k, CONST float *alpha,
                   CONST float *a, CONST INTEGER *lda, CONST float *beta,
                   float *c, CONST INTEGER *ldc);
extern VOID SSYR2K(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                   CONST INTEGER *n, CONST INTEGER *k, CONST float *alpha,
                   CONST float *a, CONST INTEGER *lda, CONST float *b,
                   CONST INTEGER *ldb, CONST float *beta, float *c,
                   CONST INTEGER *ldc);
extern VOID  STRMM(CONST CHARACTER *side, CONST CHARACTER *uplo,
                   CONST CHARACTER *transa, CONST CHARACTER *diag,
                   CONST INTEGER *m, CONST INTEGER *n, CONST float *alpha,
                   CONST float *a, CONST INTEGER *lda, float *b,
                   CONST INTEGER *ldb);
extern VOID  STRSM(CONST CHARACTER *side, CONST CHARACTER *uplo,
                   CONST CHARACTER *transa, CONST CHARACTER *diag,
                   CONST INTEGER *m, CONST INTEGER *n, CONST float *alpha,
                   CONST float *a, CONST INTEGER *lda, float *b,
                   CONST INTEGER *ldb);

#if (FORTRAN_STYLE != 1)
# define DGEMM   BLAS_NAME(dgemm,DGEMM)
# define DSYMM   BLAS_NAME(dsymm,DSYMM)
# define DSYRK   BLAS_NAME(dsyrk,DSYRK)
# define DSYR2K  BLAS_NAME(dsyr2k,DSYR2K)
# define DTRMM   BLAS_NAME(dtrmm,DTRMM)
# define DTRSM   BLAS_NAME(dtrsm,DTRSM)
#endif /* FORTRAN_STYLE */

extern VOID  DGEMM(CONST CHARACTER *transa, CONST CHARACTER *transb,
                   CONST INTEGER *m, CONST INTEGER *n, CONST INTEGER *k,
                   CONST double *alpha, CONST double *a, CONST INTEGER *lda,
                   CONST double *b, CONST INTEGER *ldb, CONST double *beta,
                   double *c, CONST INTEGER *ldc);
extern VOID  DSYMM(CONST CHARACTER *side, CONST CHARACTER *uplo,
                   CONST INTEGER *m, CONST INTEGER *n, CONST double *alpha,
                   CONST double *a, CONST INTEGER *lda, CONST double *b,
                   CONST INTEGER *ldb, CONST double *beta, double *c,
                   CONST INTEGER *ldc);
extern VOID  DSYRK(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                   CONST INTEGER *n, CONST INTEGER *k, CONST double *alpha,
                   CONST double *a, CONST INTEGER *lda, CONST double *beta,
                   double *c, CONST INTEGER *ldc);
extern VOID DSYR2K(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                   CONST INTEGER *n, CONST INTEGER *k, CONST double *alpha,
                   CONST double *a, CONST INTEGER *lda, CONST double *b,
                   CONST INTEGER *ldb, double *beta, double *c,
                   CONST INTEGER *ldc);
extern VOID  DTRMM(CONST CHARACTER *side, CONST CHARACTER *uplo,
                   CONST CHARACTER *transa, CONST CHARACTER *diag,
                   CONST INTEGER *m, CONST INTEGER *n, CONST double *alpha,
                   CONST double *a, CONST INTEGER *lda, double *b,
                   CONST INTEGER *ldb);
extern VOID  DTRSM(CONST CHARACTER *side, CONST CHARACTER *uplo,
                   CONST CHARACTER *transa, CONST CHARACTER *diag,
                   CONST INTEGER *m, CONST INTEGER *n, CONST double *alpha,
                   CONST double *a, CONST INTEGER *lda, double *b,
                   CONST INTEGER *ldb);

#if (FORTRAN_STYLE != 1)
# define CGEMM   BLAS_NAME(cgemm,CGEMM)
# define CSYMM   BLAS_NAME(csymm,CSYMM)
# define CHEMM   BLAS_NAME(chemm,CHEMM)
# define CSYRK   BLAS_NAME(csyrk,CSYRK)
# define CSYR2K  BLAS_NAME(csyr2k,CSYR2K)
# define CHERK   BLAS_NAME(cherk,CHERK)
# define CHER2K  BLAS_NAME(cher2k,CHER2K)
# define CTRMM   BLAS_NAME(ctrmm,CTRMM)
# define CTRSM   BLAS_NAME(ctrsm,CTRSM)
#endif /* FORTRAN_STYLE */

extern VOID   CGEMM(CONST CHARACTER *transa, CONST CHARACTER *transb,
                    CONST INTEGER *m, CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX8_PTR alpha, CONST COMPLEX8_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX8_PTR b,
                    CONST INTEGER *ldb, CONST COMPLEX8_PTR beta,
                    COMPLEX8_PTR c, CONST INTEGER *ldc);
extern VOID   CSYMM(CONST CHARACTER *side, CONST CHARACTER *uplo,
                    CONST INTEGER *m, CONST INTEGER *n,
                    CONST COMPLEX8_PTR alpha, CONST COMPLEX8_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX8_PTR b,
                    CONST INTEGER *ldb, CONST COMPLEX8_PTR beta,
                    COMPLEX8_PTR c, CONST INTEGER *ldc);
extern VOID   CHEMM(CONST CHARACTER *side, CONST CHARACTER *uplo,
                    CONST INTEGER *m, CONST INTEGER *n,
                    CONST COMPLEX8_PTR alpha, CONST COMPLEX8_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX8_PTR b,
                    CONST INTEGER *ldb, CONST COMPLEX8_PTR beta,
                    COMPLEX8_PTR c, CONST INTEGER *ldc);
extern VOID   CHERK(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST INTEGER *n, CONST INTEGER *k, CONST float *alpha,
                    CONST COMPLEX8_PTR a, CONST INTEGER *lda,
                    CONST float *beta, COMPLEX8_PTR c, CONST INTEGER *ldc);
extern VOID  CHER2K(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX8_PTR alpha, CONST COMPLEX8_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX8_PTR b,
                    CONST INTEGER *ldb, CONST float *beta, COMPLEX8_PTR c,
                    CONST INTEGER *ldc);
extern VOID  CHER2K(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX8_PTR alpha,
                    CONST COMPLEX8_PTR a, CONST INTEGER *lda,
                    CONST COMPLEX8_PTR b, CONST INTEGER *ldb,
                    CONST float *beta,
                    COMPLEX8_PTR c, CONST INTEGER *ldc);
extern VOID   CSYRK(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX8_PTR alpha,
                    CONST COMPLEX8_PTR a, CONST INTEGER *lda,
                    CONST COMPLEX8_PTR beta,
                    COMPLEX8_PTR c, CONST INTEGER *ldc);
extern VOID  CSYR2K(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX8_PTR alpha, CONST COMPLEX8_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX8_PTR b,
                    CONST INTEGER *ldb, CONST COMPLEX8_PTR beta, COMPLEX8_PTR c,
                    CONST INTEGER *ldc);
extern VOID   CTRMM(CONST CHARACTER *side, CONST CHARACTER *uplo,
                    CONST CHARACTER *transa, CONST CHARACTER *diag,
                    CONST INTEGER *m, CONST INTEGER *n,
                    CONST COMPLEX8_PTR alpha,
                    CONST COMPLEX8_PTR a, CONST INTEGER *lda,
                    COMPLEX8_PTR b, CONST INTEGER *ldb);
extern VOID   CTRSM(CONST CHARACTER *side, CONST CHARACTER *uplo,
                    CONST CHARACTER *transa, CONST CHARACTER *diag,
                    CONST INTEGER *m, CONST INTEGER *n,
                    CONST COMPLEX8_PTR alpha,
                    CONST COMPLEX8_PTR a, CONST INTEGER *lda,
                    COMPLEX8_PTR b, CONST INTEGER *ldb);

#if (FORTRAN_STYLE != 1)
# define ZGEMM   BLAS_NAME(zgemm,ZGEMM)
# define ZSYMM   BLAS_NAME(zsymm,ZSYMM)
# define ZHEMM   BLAS_NAME(zhemm,ZHEMM)
# define ZSYRK   BLAS_NAME(zsyrk,ZSYRK)
# define ZSYR2K  BLAS_NAME(zsyr2k,ZSYR2K)
# define ZHERK   BLAS_NAME(zherk,ZHERK)
# define ZHER2K  BLAS_NAME(zher2k,ZHER2K)
# define ZTRMM   BLAS_NAME(ztrmm,ZTRMM)
# define ZTRSM   BLAS_NAME(ztrsm,ZTRSM)
#endif /* FORTRAN_STYLE */

extern VOID   ZGEMM(CONST CHARACTER *transa, CONST CHARACTER *transb,
                    CONST INTEGER *m, CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX16_PTR b,
                    CONST INTEGER *ldb, CONST COMPLEX16_PTR beta,
                    COMPLEX16_PTR c, CONST INTEGER *ldc);
extern VOID   ZSYMM(CONST CHARACTER *side, CONST CHARACTER *uplo,
                    CONST INTEGER *m, CONST INTEGER *n,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX16_PTR b,
                    CONST INTEGER *ldb, CONST COMPLEX16_PTR beta,
                    COMPLEX16_PTR c, CONST INTEGER *ldc);
extern VOID   ZHEMM(CONST CHARACTER *side, CONST CHARACTER *uplo,
                    CONST INTEGER *m, CONST INTEGER *n,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX16_PTR b,
                    CONST INTEGER *ldb, CONST COMPLEX16_PTR beta,
                    COMPLEX16_PTR c, CONST INTEGER *ldc);
extern VOID   ZHERK(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST INTEGER *n, CONST INTEGER *k, CONST double *alpha,
                    CONST COMPLEX16_PTR a, CONST INTEGER *lda,
                    CONST double *beta, COMPLEX16_PTR c, CONST INTEGER *ldc);
extern VOID  ZHER2K(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX16_PTR b,
                    CONST INTEGER *ldb, CONST double *beta, COMPLEX16_PTR c,
                    CONST INTEGER *ldc);
extern VOID   ZSYRK(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX16_PTR beta,
                    COMPLEX16_PTR c, CONST INTEGER *ldc);
extern VOID  ZSYR2K(CONST CHARACTER *uplo, CONST CHARACTER *trans,
                    CONST INTEGER *n, CONST INTEGER *k,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR a,
                    CONST INTEGER *lda, CONST COMPLEX16_PTR b,
                    CONST INTEGER *ldb, CONST COMPLEX16_PTR beta,
                    COMPLEX16_PTR c, CONST INTEGER *ldc);
extern VOID   ZTRMM(CONST CHARACTER *side, CONST CHARACTER *uplo,
                    CONST CHARACTER *transa, CONST CHARACTER *diag,
                    CONST INTEGER *m, CONST INTEGER *n,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR a,
                    CONST INTEGER *lda, COMPLEX16_PTR b, CONST INTEGER *ldb);
extern VOID   ZTRSM(CONST CHARACTER *side, CONST CHARACTER *uplo,
                    CONST CHARACTER *transa, CONST CHARACTER *diag,
                    CONST INTEGER *m, CONST INTEGER *n,
                    CONST COMPLEX16_PTR alpha, CONST COMPLEX16_PTR a,
                    CONST INTEGER *lda, COMPLEX16_PTR b, CONST INTEGER *ldb);

_YLPK_END_DECLS

/* 
 * Local Variables:
 * mode: C
 * tab-width: 8
 * c-basic-offset: 2
 * fill-column: 78
 * coding: utf-8
 * End:
 */
