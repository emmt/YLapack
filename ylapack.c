/*
 * ylapack.c --
 *
 * Implements Yorick interface to (C/GOTO)BLAS and LAPACK.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2011-2012 Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
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

#include <yapi.h>
#include <ydata.h>
#include <pstdlib.h>
#include <stdio.h>
#include <string.h>

#include "ylapack.h"

#define MAX(a, b)       ((a) >= (b) ? (a) : (b))
#define MIN(a, b)       ((a) <= (b) ? (a) : (b))
#define TRUE            (1)
#define FALSE           (0)

#define INCR_ADDRESS(ptr, nbytes)  \
  ((void *)((unsigned char *)(ptr) + (nbytes)))

#define ROUND_UP(N, S) ((((S) - 1 + (N))/(S))*(S))

#define YSCRATCH_BROKEN 0
#define YTEMP_BROKEN    0
#define YSWAP_BROKEN    0

#if YSCRATCH_BROKEN
static y_userobj_t worskspace_class = {
  /* type_name  */ "worskspace",
  /* on_free    */ NULL,
  /* on_print   */ NULL,
  /* on_eval    */ NULL,
  /* on_extract */ NULL,
  /* uo_ops     */ NULL
};
#define ypush_scratch(size, unused) ypush_obj(&worskspace_class, size)
#endif

static long get_size_of_linear_system(const long adims[], const long bdims[]);
/* Checks that left hand side matrix A and right hand side vector B
   have compatible dimensions.  Returns the number of linear equations. */

static float  get_f(int iarg);
static double get_d(int iarg);
static void   get_z(int iarg, double z[2]);
/* Fetch a float/double/complex scalar form Yorick stack. */

static void push_string(const char *value);

typedef void *push_array_t(long *dims);
static  float *push_copy_f(const  float *src, const long dims[]);
static double *push_copy_d(const double *src, const long dims[]);
static double *push_copy_z(const double *src, const long dims[]);

static void get_lapack_trans(int iarg, CHARACTER trans[2]);
static void get_lapack_uplo(int iarg, CHARACTER uplo[2]);

/* The following structure is used to collect all the informations
   representing a linear system like A.X = B. */
typedef struct _linear_system linear_system_t;
struct _linear_system {
  void *a, *b; /* pointers to LHS matrix and RHS vector */
  int type; /* data type of A and B (Y_DOUBLE, Y_FLOAT or Y_COMPLEX) */
  INTEGER n, nrhs, lda, ldb; /* system dimensions */
};

static long get_linear_system(linear_system_t *sys,
                              int a_iarg, int b_iarg,
                              unsigned int flags);
/* Extracts the components of a linear system of equations from Yorick stack.
   A_IARG and B_IARG are the stack indices of the LHS matrix A and of the RHS
   vector B, FLAGS is used to check the type of the arguments (and possibly to
   convert them) and is one of: */
#define REAL_SYSTEM            (1 << 0)
#define COMPLEX_SYSTEM         (1 << 1)
#define REAL_OR_COMPLEX_SYSTEM (REAL_SYSTEM | COMPLEX_SYSTEM)


#ifdef USE_CBLAS
static CBLAS_DIAG get_cblas_diag(int iarg);
static CBLAS_ORDER get_cblas_order(int iarg);
static CBLAS_SIDE get_cblas_side(int iarg);
static CBLAS_TRANSPOSE get_cblas_trans(int iarg);
static CBLAS_UPLO get_cblas_uplo(int iarg);
#endif

#ifdef USE_CBLAS
# define BLAS_ALTERNATIVE(cblas_expr, blas_expr) (cblas_expr)
# define GET_BLAS_DIAG(iarg, var)   var = get_cblas_diag(iarg)
# define GET_BLAS_ORDER(iarg, var)  var = get_cblas_order(iarg)
# define GET_BLAS_SIDE(iarg, var)   var = get_cblas_side(iarg)
# define GET_BLAS_TRANS(iarg, var)  var = get_cblas_trans(iarg)
# define GET_BLAS_UPLO(iarg, var)   var = get_cblas_uplo(iarg)
#else
# define BLAS_ALTERNATIVE(cblas_expr, blas_expr) (blas_expr)
# define GET_BLAS_DIAG(iarg, var)   get_lapack_diag(iarg, var)
# define GET_BLAS_ORDER(iarg, var)  get_lapack_order(iarg, var)
# define GET_BLAS_SIDE(iarg, var)   get_lapack_side(iarg, var)
# define GET_BLAS_TRANS(iarg, var)  get_lapack_trans(iarg, var)
# define GET_BLAS_UPLO(iarg, var)   get_lapack_uplo(iarg, var)
#endif


/* ORDER: */
#define LPK_ROW_MAJOR  101
#define LPK_COL_MAJOR  102
/* TRANSPOSE: */
#define LPK_NO_TRANS   111
#define LPK_TRANS      112
#define LPK_CONJ_TRANS 113
#define LPK_ATLAS_CONJ 114
/* UPLO: */
#define LPK_UPPER      121
#define LPK_LOWER      122
/* DIAG: */
#define LPK_NON_UNIT   131
#define LPK_UNIT       132
/* SIDE: */
#define LPK_LEFT       141
#define LPK_RIGHT      142

/* Make some basic checking. (Needed for tables of functions, etc.) */
#if (Y_CHAR != 0)
# error assertion failed (Y_CHAR != 0)
#endif
#if (Y_SHORT != 1)
# error assertion failed (Y_SHORT != 1)
#endif
#if (Y_INT != 2)
# error assertion failed (Y_INT != 2)
#endif
#if (Y_LONG != 3)
# error assertion failed (Y_LONG != 3)
#endif
#if (Y_FLOAT != 4)
# error assertion failed (Y_FLOAT != 4)
#endif
#if (Y_DOUBLE != 5)
# error assertion failed (Y_DOUBLE != 5)
#endif
#if (Y_COMPLEX != 6)
# error assertion failed (Y_COMPLEX != 6)
#endif

/* PRIVATE DATA */
static long full_index = -1L;
static long slow_index = -1L;
static char msgbuf[128]; /* used for error messages */

/*---------------------------------------------------------------------------*/
/* HELPER MACROS */

/* There are many repeated tasks in parsing the arguments from Yorick
   stack.  Some of these tasks are implemented by helper macros.

   The following helper macros are used to fetch informations from Yorick
   stack.  For a symbol "sym", the following variables are used:
   void *sym_ptr;  // address of array
   long sym_ntot;  // total number of elements
   long sym_cnt;   // number of elements to consider
   long sym_off;   // 0-based offset to 1st element to consider
   long sym_stp;   // increment between elements of interest
   int sym_type;   // data type (e.g. Y_DOUBLE)
   int sym_iarg;   // stack index of symbol
*/

/* The macro GET_ARRAY gets an array for symbol "sym", NAME is for error
   message and MAX_TYPE is Y_DOUBLE or Y_COMPLEX. */
#define GET_ARRAY(sym, name, max_type) do {                     \
    sym##_ptr = ygeta_any(sym##_iarg, &sym##_ntot, sym##_dims,  \
                          &sym##_type);                         \
    if (sym##_type < Y_CHAR || sym##_type > max_type) {         \
      y_error("bad data type for " name);                       \
    }                                                           \
  } while (0)

/* The macro GET_RANGE parses the range at position IARG and update variables
   sym_ntot, sym_step and sym_offset which must exists and sym_ntot must have
   been set with the number of element in symbol. */
#define GET_RANGE(sym, iarg) do {               \
    int temp_iarg = (iarg);                     \
    if (temp_iarg >= 0) {                       \
      long temp_range[3];                       \
      parse_range(RANGE_CBLAS_STYLE, temp_iarg, \
                  sym##_ntot, temp_range);      \
      sym##_off = temp_range[0];                \
      sym##_stp = temp_range[1];                \
      sym##_cnt = temp_range[2];                \
    } else {                                    \
      sym##_off = 0;                            \
      sym##_stp = 1;                            \
      sym##_cnt = sym##_ntot;                   \
    }                                           \
  } while (0)

/* Convert data type. */
#define COERCE(sym, type)   do {                                \
    sym##_ptr = ygeta_coerce(sym##_iarg, sym##_ptr, sym##_ntot, \
                             sym##_dims, sym##_type, type);     \
    sym##_type = type;                                          \
  } while (0)

/*---------------------------------------------------------------------------*/
/* BUILT-IN UTILITIES */

void Y_lpk_init(int argc)
{
  ypush_long(1);
  ypush_long(2);
  yarg_swap(0, 1);
  if (ygets_l(0) != 1 || ygets_l(1) != 2) {
    y_error("yarg_swap broken in yapi.c");
  }
  full_index = yget_global("full", 0);
  slow_index = yget_global("slow", 0);
  ypush_nil();
}

void Y_lpk_model(int argc)
{
#if defined(USE_MKL)
  push_string("MKL (Intel Math kernel library)");
#elif defined(USE_GOTOBLAS)
  push_string("GotoBLAS2");
#elif defined(USE_ATLAS)
  push_string("Atlas");
#else
  push_string("Lapack");
#endif
}

void Y_lpk_laver(int argc)
{
  long dims[2], *result;
  INTEGER major, minor, patch;

  ILAVER(&major, &minor, &patch);
  dims[0] = 1;
  dims[1] = 3;
  result = ypush_l(dims);
  result[0] = major;
  result[1] = minor;
  result[2] = patch;
}

/*---------------------------------------------------------------------------*/
/* LEVEL 1 BLAS */

#ifdef USE_CBLAS
# define BLAS_INTEGER CBLAS_INTEGER
#else
# define BLAS_INTEGER INTEGER
#endif

typedef struct  _vector  vector_t;
typedef struct  _matrix  matrix_t;
typedef struct  _array   array_t;
typedef struct  _scalar  scalar_t;

static long elem_size[] = {
  sizeof(char),
  sizeof(short),
  sizeof(int),
  sizeof(long),
  sizeof(float),
  sizeof(double),
  sizeof(double)*2
};

struct _scalar {
  union {
    char c;
    short s;
    int i;
    long l;
    float f;
    double d;
    double z[2];
  } value;
  int type;
};

/* Structure to store a general Yorick array. */
struct _array {
  void *data;           /* address of first element */
  long ntot;            /* total number of elements */
  long dims[Y_DIMSIZE]; /* dimension list, dims[0] is the rank */
  int type;             /* data type identifier */
  int iarg;             /* index in argument list; used for conversion
                           (*beware* of stack changes when
                           pushing/dropping/swapping arguments) */
};

/* A "vector" is a wrapper around a Yorick array intended to represent a
   vector for BLAS calls.  It may start at another element than the first one
   and have a non-unit increment. */
struct _vector {
  array_t arr;       /* Yorick array */
  long offset;       /* 0-based index of first element to consider */
  BLAS_INTEGER size; /* number of elements to consider */
  BLAS_INTEGER step; /* increment between successive elements to consider (can
                        be negative) */
  int rank;          /* rank of the "vector", 1 for a flat vector; otherwise
                        can be greater than 1 */
};

/* Macros for quick access to vector contents. (*Beware* that argument may be
   used more than once.) */
#define VECTOR_RANK(vec) ((vec).rank)
#define VECTOR_SIZE(vec) ((vec).size)
#define VECTOR_STEP(vec) ((vec).step)
#define VECTOR_TYPE(vec) ((vec).arr.type)
#define VECTOR_BASE(vec) INCR_ADDRESS((vec).arr.data, \
                                   elem_size[(vec).arr.type]*(vec).offset)

/* A "matrix" is a wrapper around a Yorick array intended to represent a
   matrix for BLAS or LAPACK calls. */
struct _matrix {
  array_t arr;          /* Yorick array */
  INTEGER nrows, ncols; /* number of rows (leading dimension) and columns */
};

/* Macros for quick access to matrix contents.  (The leading dimension of a
   matrix is its number of rows.) */
#define MATRIX_TYPE(mat)   ((mat).arr.type)
#define MATRIX_BASE(mat)   ((mat).arr.data)
#define MATRIX_RANK(mat)   ((mat).arr.dims[0])
#define MATRIX_NROWS(mat)  ((mat).nrows)
#define MATRIX_NCOLS(mat)  ((mat).ncols)
#define MATRIX_STRIDE(mat) MATRIX_NROWS(mat)

/* Functions to fetch array-based objects from Yorick stack.  IARG is the
   stack index of the object; for sub-indexing vectors, RNG is the stack index
   of the range, -1 if none. */
static void get_array(array_t *arr, int iarg);
static void get_vector(vector_t *vec, int iarg, int rng);
static void get_matrix(matrix_t *mat, int iarg);

/* Functions to convert array-based objects to a given data type (*beware* of
   stack changes when pushing/dropping/swapping arguments). */
static void coerce_array(array_t *arr, int type);
static void coerce_vector(vector_t *vec, int type);
static void coerce_matrix(matrix_t *mat, int type);

static BLAS_INTEGER check_vector_sizes(const vector_t *x, const vector_t *y);

static void copy_array(array_t *arr);
/* Replace stack argument arr->iarg by a copy of the array. */

/* The BLAS_ARG_* and LAPACK_ARG_* macros can be used to put an argument in
 * the argument list for a call to a BLAS and LAPACK function respectively.
 *   BLAS_ARG_I / LAPACK_ARG_I - for an integer scalar argument;
 *   BLAS_ARG_R / LAPACK_ARG_R - for a real scalar argument;
 *   BLAS_ARG_C / LAPACK_ARG_C - for a complex scalar argument (a pair of
 *                               reals);
 *   BLAS_ARG_V                - for a BLAS vector argument (array address
 *                               and increment);
 *   BLAS_ARG_M / LAPACK_ARG_M - for a matrix argument (array address and
 *                               size of leading dimension);
 */
#ifdef USE_CBLAS
# define BLAS_ARG_I(i)     (i)
# define BLAS_ARG_R(r)     (r)
# define BLAS_ARG_C(c)     (c)
#else
# define BLAS_ARG_I(i)    &(i)
# define BLAS_ARG_R(r)    &(r)
# define BLAS_ARG_C(c)     (c)
#endif
#define LAPACK_ARG_I(i)   &(i)
#define LAPACK_ARG_R(r)   &(r)
#define LAPACK_ARG_C(c)    (c)
#define BLAS_ARG_V(vec)   VECTOR_BASE(vec), BLAS_ARG_I(VECTOR_STEP(vec))
#define BLAS_ARG_M(mat)   MATRIX_BASE(mat), BLAS_ARG_I(MATRIX_STRIDE(mat))
#define LAPACK_ARG_M(mat) MATRIX_BASE(mat), LAPACK_ARG_I(MATRIX_STRIDE(mat))


/* Macros for BLAS level 1

   CALL_...(fn,cn,...)
   FN is the FORTRAN name of the function to call
   CN is the C name of the function to call
   The name is built from the type of arguments:
   I = integer argument, F = flag,
   R = real scalar, C = complex scalar,
   V = vector, M = matrix, A = array
 */

#ifdef USE_CBLAS
# define CALL_FUNC(fn,cn) cn
#else
# define CALL_FUNC(fn,cn) fn
#endif

#define CALL_BLAS1_V(fn,cn,n,x)                 \
  CALL_FUNC(fn,cn)(BLAS_ARG_I(n),               \
                   BLAS_ARG_V(x))
#define CALL_BLAS1_RV(fn,cn,n,alpha,x)          \
  CALL_FUNC(fn,cn)(BLAS_ARG_I(n),               \
                   BLAS_ARG_R(alpha),           \
                   BLAS_ARG_V(x))
#define CALL_BLAS1_CV(fn,cn,n,alpha,x)          \
  CALL_FUNC(fn,cn)(BLAS_ARG_I(n),               \
                   BLAS_ARG_C(alpha),           \
                   BLAS_ARG_V(x))
#define CALL_BLAS1_VV(fn,cn,n,x,y)              \
  CALL_FUNC(fn,cn)(BLAS_ARG_I(n),               \
                   BLAS_ARG_V(x),               \
                   BLAS_ARG_V(y))
#define CALL_BLAS1_RVV(fn,cn,n,alpha,x,y)       \
  CALL_FUNC(fn,cn)(BLAS_ARG_I(n),               \
                   BLAS_ARG_R(alpha),           \
                   BLAS_ARG_V(x),               \
                   BLAS_ARG_V(y))
#define CALL_BLAS1_CVV(fn,cn,n,alpha,x,y)       \
  CALL_FUNC(fn,cn)(BLAS_ARG_I(n),               \
                   BLAS_ARG_C(alpha),           \
                   BLAS_ARG_V(x),               \
                   BLAS_ARG_V(y))

#define CALL_SSCAL(n,alpha,x) CALL_BLAS1_RV(SSCAL,cblas_sscal,n,alpha,x)
#define CALL_DSCAL(n,alpha,x) CALL_BLAS1_RV(DSCAL,cblas_dscal,n,alpha,x)
#define CALL_CSCAL(n,alpha,x) CALL_BLAS1_CV(CSCAL,cblas_cscal,n,alpha,x)
#define CALL_ZSCAL(n,alpha,x) CALL_BLAS1_CV(ZSCAL,cblas_zscal,n,alpha,x)

#define CALL_SSWAP(n,x,y) CALL_BLAS1_VV(SSWAP,cblas_sswap,n,x,y)
#define CALL_DSWAP(n,x,y) CALL_BLAS1_VV(DSWAP,cblas_dswap,n,x,y)
#define CALL_CSWAP(n,x,y) CALL_BLAS1_VV(CSWAP,cblas_cswap,n,x,y)
#define CALL_ZSWAP(n,x,y) CALL_BLAS1_VV(ZSWAP,cblas_zswap,n,x,y)

#define CALL_SCOPY(n,x,y) CALL_BLAS1_VV(SCOPY,cblas_scopy,n,x,y)
#define CALL_DCOPY(n,x,y) CALL_BLAS1_VV(DCOPY,cblas_dcopy,n,x,y)
#define CALL_CCOPY(n,x,y) CALL_BLAS1_VV(CCOPY,cblas_ccopy,n,x,y)
#define CALL_ZCOPY(n,x,y) CALL_BLAS1_VV(ZCOPY,cblas_zcopy,n,x,y)

#define CALL_SAXPY(n,alpha,x,y) CALL_BLAS1_RVV(SAXPY,cblas_saxpy,n,alpha,x,y)
#define CALL_DAXPY(n,alpha,x,y) CALL_BLAS1_RVV(DAXPY,cblas_daxpy,n,alpha,x,y)
#define CALL_CAXPY(n,alpha,x,y) CALL_BLAS1_CVV(CAXPY,cblas_caxpy,n,alpha,x,y)
#define CALL_ZAXPY(n,alpha,x,y) CALL_BLAS1_CVV(ZAXPY,cblas_zaxpy,n,alpha,x,y)

#define CALL_SASUM(n,x)  CALL_BLAS1_V(SASUM,cblas_sasum,n,x)
#define CALL_DASUM(n,x)  CALL_BLAS1_V(DASUM,cblas_dasum,n,x)
#define CALL_SCASUM(n,x) CALL_BLAS1_V(SCASUM,cblas_scasum,n,x)
#define CALL_DZASUM(n,x) CALL_BLAS1_V(DZASUM,cblas_dzasum,n,x)

#define CALL_SNRM2(n,x)  CALL_BLAS1_V(SNRM2,cblas_snrm2,n,x)
#define CALL_DNRM2(n,x)  CALL_BLAS1_V(DNRM2,cblas_dnrm2,n,x)
#define CALL_SCNRM2(n,x) CALL_BLAS1_V(SCNRM2,cblas_scnrm2,n,x)
#define CALL_DZNRM2(n,x) CALL_BLAS1_V(DZNRM2,cblas_dznrm2,n,x)

#ifdef USE_CBLAS
# define CALL_DOT_C(fn,cn,r,n,x,y)              \
  cn(BLAS_ARG_I(n),                             \
     BLAS_ARG_V(x),                             \
     BLAS_ARG_V(y),                             \
     (void *)r)
#else
# define CALL_DOT_C(fn,cn,r,n,x,y)              \
  fn((void *)r,                                 \
     BLAS_ARG_I(n),                             \
     BLAS_ARG_V(x),                             \
     BLAS_ARG_V(y))
#endif
#define CALL_SDOT(n,x,y)    CALL_BLAS1_VV(SDOT,cblas_sdot,n,x,y)
#define CALL_DDOT(n,x,y)    CALL_BLAS1_VV(DDOT,cblas_ddot,n,x,y)
#define CALL_CDOTU(r,n,x,y) CALL_DOT_C(CDOTU, cblas_cdotu_sub, r, n, x, y)
#define CALL_CDOTC(r,n,x,y) CALL_DOT_C(CDOTC, cblas_cdotc_sub, r, n, x, y)
#define CALL_ZDOTU(r,n,x,y) CALL_DOT_C(ZDOTU, cblas_zdotu_sub, r, n, x, y)
#define CALL_ZDOTC(r,n,x,y) CALL_DOT_C(ZDOTC, cblas_zdotc_sub, r, n, x, y)

#ifdef USE_CBLAS
# define CALL_GEMV_R(fn,cn,trans,m,n,alpha,A,x,beta,y)  \
  cn(ORDER, trans, BLAS_ARG_I(m), BLAS_ARG_I(n),        \
     BLAS_ARG_R(alpha),  BLAS_ARG_M(A), BLAS_ARG_V(x),  \
     BLAS_ARG_R(beta), BLAS_ARG_V(y))
# define CALL_GEMV_C(fn,cn,trans,m,n,alpha,A,x,beta,y)  \
  cn(ORDER, trans, BLAS_ARG_I(m), BLAS_ARG_I(n),        \
     BLAS_ARG_C(alpha),  BLAS_ARG_M(A), BLAS_ARG_V(x),  \
     BLAS_ARG_C(beta), BLAS_ARG_V(y))
# define CALL_SYMV_R(fn,cn,uplo,n,alpha,A,x,beta,y)     \
  cn(ORDER, uplo, BLAS_ARG_I(n),                        \
     BLAS_ARG_R(alpha),  BLAS_ARG_M(A), BLAS_ARG_V(x),  \
     BLAS_ARG_R(beta), BLAS_ARG_V(y))
# define CALL_HEMV_C(fn,cn,uplo,n,alpha,A,x,beta,y)     \
  cn(ORDER, uplo, BLAS_ARG_I(n),                        \
     BLAS_ARG_C(alpha),  BLAS_ARG_M(A), BLAS_ARG_V(x),  \
     BLAS_ARG_C(beta), BLAS_ARG_V(y))
#else
# define CALL_GEMV_R(fn,cn,trans,m,n,alpha,A,x,beta,y)  \
  fn(trans, BLAS_ARG_I(m), BLAS_ARG_I(n),               \
     BLAS_ARG_R(alpha),  BLAS_ARG_M(A), BLAS_ARG_V(x),  \
     BLAS_ARG_R(beta), BLAS_ARG_V(y))
# define CALL_GEMV_C(fn,cn,trans,m,n,alpha,A,x,beta,y)  \
  fn(trans, BLAS_ARG_I(m), BLAS_ARG_I(n),               \
     BLAS_ARG_C(alpha),  BLAS_ARG_M(A), BLAS_ARG_V(x),  \
     BLAS_ARG_C(beta), BLAS_ARG_V(y))
# define CALL_SYMV_R(fn,cn,uplo,n,alpha,A,x,beta,y)     \
  fn(uplo, BLAS_ARG_I(n),                               \
     BLAS_ARG_R(alpha),  BLAS_ARG_M(A), BLAS_ARG_V(x),  \
     BLAS_ARG_R(beta), BLAS_ARG_V(y))
# define CALL_HEMV_C(fn,cn,uplo,n,alpha,A,x,beta,y)     \
  fn(uplo, BLAS_ARG_I(n),                               \
     BLAS_ARG_C(alpha),  BLAS_ARG_M(A), BLAS_ARG_V(x),  \
     BLAS_ARG_C(beta), BLAS_ARG_V(y))
#endif

#define CALL_SGEMV(trans,m,n,alpha,A,x,beta,y) \
  CALL_GEMV_R(SGEMV,cblas_sgemv,trans,m,n,alpha,A,x,beta,y)
#define CALL_DGEMV(trans,m,n,alpha,A,x,beta,y) \
  CALL_GEMV_R(DGEMV,cblas_dgemv,trans,m,n,alpha,A,x,beta,y)

#define CALL_CGEMV(trans,m,n,alpha,A,x,beta,y) \
  CALL_GEMV_C(CGEMV,cblas_cgemv,trans,m,n,alpha,A,x,beta,y)
#define CALL_ZGEMV(trans,m,n,alpha,A,x,beta,y) \
  CALL_GEMV_C(ZGEMV,cblas_zgemv,trans,m,n,alpha,A,x,beta,y)

#define CALL_SSYMV(uplo,n,alpha,A,x,beta,y) \
  CALL_SYMV_R(SSYMV,cblas_ssymv,uplo,n,alpha,A,x,beta,y)
#define CALL_DSYMV(uplo,n,alpha,A,x,beta,y) \
  CALL_SYMV_R(DSYMV,cblas_dsymv,uplo,n,alpha,A,x,beta,y)

#define CALL_CHEMV(uplo,n,alpha,A,x,beta,y) \
  CALL_HEMV_C(CHEMV,cblas_chemv,uplo,n,alpha,A,x,beta,y)
#define CALL_ZHEMV(uplo,n,alpha,A,x,beta,y) \
  CALL_HEMV_C(ZHEMV,cblas_zhemv,uplo,n,alpha,A,x,beta,y)

#ifdef USE_CBLAS
# define CALL_GEMM_R(fn,cn,transa,transb,m,n,k,alpha,A,B,beta,C)        \
  cn(ORDER,transa,transb,BLAS_ARG_I(m),BLAS_ARG_I(n),BLAS_ARG_I(k),     \
     BLAS_ARG_R(alpha),BLAS_ARG_M(A),BLAS_ARG_M(B),                     \
     BLAS_ARG_R(beta),BLAS_ARG_M(C))
# define CALL_GEMM_C(fn,cn,transa,transb,m,n,k,alpha,A,B,beta,C)        \
  cn(ORDER,transa,transb,BLAS_ARG_I(m),BLAS_ARG_I(n),BLAS_ARG_I(k),     \
     BLAS_ARG_C(alpha),BLAS_ARG_M(A),BLAS_ARG_M(B),                     \
     BLAS_ARG_C(beta),BLAS_ARG_M(C))
#else
# define CALL_GEMM_R(fn,cn,transa,transb,m,n,k,alpha,A,B,beta,C)        \
  fn(transa,transb,BLAS_ARG_I(m),BLAS_ARG_I(n),BLAS_ARG_I(k),           \
     BLAS_ARG_R(alpha),BLAS_ARG_M(A),BLAS_ARG_M(B),                     \
     BLAS_ARG_R(beta),BLAS_ARG_M(C))
# define CALL_GEMM_C(fn,cn,transa,transb,m,n,k,alpha,A,B,beta,C)        \
  fn(transa,transb,BLAS_ARG_I(m),BLAS_ARG_I(n),BLAS_ARG_I(k),           \
     BLAS_ARG_C(alpha),BLAS_ARG_M(A),BLAS_ARG_M(B),                     \
     BLAS_ARG_C(beta),BLAS_ARG_M(C))
#endif

#define CALL_SGEMM(transa,transb,m,n,k,alpha,A,B,beta,C)                \
  CALL_GEMM_R(SGEMM,cblas_sgemm,transa,transb,m,n,k,alpha,A,B,beta,C)
#define CALL_DGEMM(transa,transb,m,n,k,alpha,A,B,beta,C)                \
  CALL_GEMM_R(DGEMM,cblas_dgemm,transa,transb,m,n,k,alpha,A,B,beta,C)
#define CALL_CGEMM(transa,transb,m,n,k,alpha,A,B,beta,C)                \
  CALL_GEMM_C(CGEMM,cblas_cgemm,transa,transb,m,n,k,alpha,A,B,beta,C)
#define CALL_ZGEMM(transa,transb,m,n,k,alpha,A,B,beta,C)                \
  CALL_GEMM_C(ZGEMM,cblas_zgemm,transa,transb,m,n,k,alpha,A,B,beta,C)

#ifdef USE_CBLAS
# define _CALL_SYRK(fn,cn,uplo,trans,n,k,alpha,A,beta,C)                \
  cn(ORDER,uplo,trans,BLAS_ARG_I(n),BLAS_ARG_I(k),BLAS_ARG_R(alpha),    \
     BLAS_ARG_M(A),BLAS_ARG_R(beta),BLAS_ARG_M(C))
#else
# define _CALL_SYRK(fn,cn,uplo,trans,n,k,alpha,A,beta,C)                \
  fn(uplo,trans,BLAS_ARG_I(n),BLAS_ARG_I(k),BLAS_ARG_R(alpha),          \
     BLAS_ARG_M(A),BLAS_ARG_R(beta),BLAS_ARG_M(C))
#endif
#define CALL_SSYRK(uplo,trans,n,k,alpha,A,beta,C)               \
  _CALL_SYRK(SSYRK,cblas_ssyrk,uplo,trans,n,k,alpha,A,beta,C)
#define CALL_DSYRK(uplo,trans,n,k,alpha,A,beta,C)               \
  _CALL_SYRK(DSYRK,cblas_dsyrk,uplo,trans,n,k,alpha,A,beta,C)

void Y_lpk_nrm2(int argc)
{
  vector_t x;
  BLAS_INTEGER n;

  if (argc < 1 || argc > 2) {
    y_error("lpk_nrm2 takes 1 or 2 arguments");
  }
  get_vector(&x, argc - 1, argc - 2);
  if (VECTOR_TYPE(x) < Y_FLOAT) {
    coerce_vector(&x, Y_FLOAT);
  }
  n = x.size;
  if (VECTOR_TYPE(x) == Y_DOUBLE) {
    ypush_double(CALL_DNRM2(n,x));
  } else if (VECTOR_TYPE(x) == Y_FLOAT) {
    ypush_f(NULL)[0] = CALL_SNRM2(n,x);
  } else {
    ypush_double(CALL_DZNRM2(n,x));
  }
}

void Y_lpk_asum(int argc)
{
  vector_t x;
  BLAS_INTEGER n;

  if (argc < 1 || argc > 2) {
    y_error("lpk_asum takes 1 or 2 arguments");
  }
  get_vector(&x, argc - 1, argc - 2);
  if (VECTOR_TYPE(x) < Y_FLOAT) {
    coerce_vector(&x, Y_FLOAT);
  }
  n = x.size;
  if (VECTOR_TYPE(x) == Y_DOUBLE) {
    ypush_double(CALL_DASUM(n,x));
  } else if (VECTOR_TYPE(x) == Y_FLOAT) {
    ypush_f(NULL)[0] = CALL_SASUM(n,x);
  } else {
    ypush_double(CALL_DZASUM(n,x));
  }
}

void Y_lpk_copy(int argc)
{
  vector_t x, y;
  long y_ref;
  BLAS_INTEGER n;
  int drop, type;

  if (argc < 2 || argc > 4) {
    y_error("lpk_copy takes from 2 to 4 arguments");
  }
  if (argc == 2) {
    y_ref = yget_ref(argc - 2);
    get_vector(&x, argc - 1, -1);
    get_vector(&y, argc - 2, -1);
    drop = FALSE;
  } else if (argc == 3) {
    y_ref = yget_ref(argc - 2); /* get ref *before* yarg_typeid */
    if (yarg_typeid(argc - 2) == Y_RANGE) {
      y_ref = yget_ref(argc - 3);
      get_vector(&x, argc - 1, argc - 2);
      get_vector(&y, argc - 3, -1);
      drop = FALSE;
    } else {
      get_vector(&x, argc - 1, -1);
      get_vector(&y, argc - 2, argc - 3);
      drop = TRUE;
    }
  } else {
    y_ref = yget_ref(argc - 3);
    get_vector(&x, argc - 1, argc - 2);
    get_vector(&y, argc - 3, argc - 4);
    drop = TRUE;
  }
  if (y_ref < 0) {
    y_error("argument Y must be a variable");
  }
  n = check_vector_sizes(&x, &y);
  type = Y_FLOAT;
  if (type < VECTOR_TYPE(x)) {
    type = VECTOR_TYPE(x);
  }
  if (type < VECTOR_TYPE(y)) {
    type = VECTOR_TYPE(y);
  }
  if (VECTOR_TYPE(x) != type) {
    coerce_vector(&x, type);
  }
  if (VECTOR_TYPE(y) != type) {
    coerce_vector(&y, type);
  } else {
    y_ref = -1L; /* no needs to re-define global symbol */
  }
  if (type == Y_DOUBLE) {
    CALL_DCOPY(n,x,y);
  } else if (type == Y_FLOAT) {
    CALL_SCOPY(n,x,y);
  } else {
    CALL_ZCOPY(n,x,y);
  }

  /* Update global variable and make sure Y is left on top of the stack. */
  if (y_ref >= 0L) {
    yput_global(y_ref, y.arr.iarg);
  }
  if (drop) {
    yarg_drop(1);
  }
}

void Y_lpk_swap(int argc)
{
  vector_t x, y;
  long x_ref, y_ref;
  BLAS_INTEGER n;
  int type;

  if (argc < 2 || argc > 4) {
    y_error("lpk_swap takes from 2 to 4 arguments");
  }
  x_ref = yget_ref(argc - 1);
  if (x_ref < 0) {
    y_error("argument X must be a variable");
  }
  if (argc == 2) {
    y_ref = yget_ref(argc - 2);
    get_vector(&x, argc - 1, -1);
    get_vector(&y, argc - 2, -1);
  } else if (argc == 3) {
    y_ref = yget_ref(argc - 2); /* get ref *before* yarg_typeid */
    if (yarg_typeid(argc - 2) == Y_RANGE) {
      y_ref = yget_ref(argc - 3);
      get_vector(&x, argc - 1, argc - 2);
      get_vector(&y, argc - 3, -1);
    } else {
      get_vector(&x, argc - 1, -1);
      get_vector(&y, argc - 2, argc - 3);
    }
  } else {
    y_ref = yget_ref(argc - 3);
    get_vector(&x, argc - 1, argc - 2);
    get_vector(&y, argc - 3, argc - 4);
  }
  if (y_ref < 0) {
    y_error("argument Y must be a variable");
  }
  n = check_vector_sizes(&x, &y);
  type = Y_FLOAT;
  if (type < VECTOR_TYPE(x)) {
    type = VECTOR_TYPE(x);
  }
  if (type < VECTOR_TYPE(y)) {
    type = VECTOR_TYPE(y);
  }
  if (VECTOR_TYPE(x) != type) {
    coerce_vector(&x, type);
  } else {
    x_ref = -1L; /* no needs to re-define global symbol */
  }
  if (VECTOR_TYPE(y) != type) {
    coerce_vector(&y, type);
  } else {
    y_ref = -1L; /* no needs to re-define global symbol */
  }
  if (type == Y_DOUBLE) {
    CALL_DSWAP(n,x,y);
  } else if (type == Y_FLOAT) {
    CALL_SSWAP(n,x,y);
  } else {
    CALL_ZSWAP(n,x,y);
  }
  if (x_ref >= 0L) {
    yput_global(x_ref, x.arr.iarg);
  }
  if (y_ref >= 0L) {
    yput_global(y_ref, y.arr.iarg);
  }
  ypush_nil();
}

void Y_lpk_scal(int argc)
{
  vector_t x;
  long x_ref;
  BLAS_INTEGER n;
  int type;

  if (argc < 2 || argc > 3) {
    y_error("lpk_scal takes 2 or 3 arguments");
  }
  x_ref = yget_ref(argc - 2);
  if (x_ref < 0) {
    y_error("argument X must be a variable");
  }
  get_vector(&x, argc - 2, argc - 3);
  if (VECTOR_TYPE(x) < Y_FLOAT) {
    coerce_vector(&x, Y_FLOAT);
  } else {
    x_ref = -1L; /* no needs to re-define global symbol */
  }
  type = VECTOR_TYPE(x);
  n = x.size;
  if (type == Y_DOUBLE) {
    double alpha = get_d(argc - 1);
    CALL_DSCAL(n,alpha,x);
  } else if (type == Y_FLOAT) {
    float alpha = get_f(argc - 1);
    CALL_SSCAL(n,alpha,x);
  } else {
    double alpha[2];
    get_z(argc - 1, alpha);
    CALL_ZSCAL(n,alpha,x);
  }

  /* If X has been converted, save into its variable and manage to left X on
     top of the satck. */
  if (x_ref >= 0L) {
    yput_global(x_ref, x.arr.iarg);
  }
  if (argc > 2) {
    yarg_drop(1);
  }
}

void Y_lpk_dot(int argc)
{
  vector_t x, y;
  int type;
  BLAS_INTEGER n;

  /* Get the arguments. */
  if (argc == 2) {
    get_vector(&x, argc - 1, -1);
    get_vector(&y, argc - 2, -1);
  } else if (argc == 3) {
    if (yarg_typeid(argc - 2) == Y_RANGE) {
      get_vector(&x, argc - 1, argc - 2);
      get_vector(&y, argc - 3, -1);
    } else {
      get_vector(&x, argc - 1, -1);
      get_vector(&y, argc - 2, argc - 3);
    }
  } else if (argc == 4) {
    get_vector(&x, argc - 1, argc - 2);
    get_vector(&y, argc - 3, argc - 4);
 } else {
    y_error("expecting 2, 3 or 4 arguments");
    return;
  }

  /* Check sizes and types. */
  n = check_vector_sizes(&x, &y);
  type = Y_FLOAT;
  if (type < VECTOR_TYPE(x)) {
    type = VECTOR_TYPE(x);
  }
  if (type < VECTOR_TYPE(y)) {
    type = VECTOR_TYPE(y);
  }
  if (VECTOR_TYPE(x) != type) {
    coerce_vector(&x, type);
  }
  if (VECTOR_TYPE(y) != type) {
    coerce_vector(&y, type);
  }

  /* Apply the operation. */
  if (type == Y_DOUBLE) {
    ypush_double(CALL_DDOT(n,x,y));
  } else if (type == Y_FLOAT) {
    *ypush_f(NULL) = CALL_DDOT(n,x,y);
  } else {
    double *result = ypush_z(NULL);
    CALL_ZDOTC(result, n, x, y);
  }
}

void Y_lpk_axpy(int argc)
{
  vector_t x, y;
  long y_ref;
  BLAS_INTEGER n;
  int type, drop;

  /* Get the vector arguments. */
  if (argc == 3) {
    y_ref = yget_ref(argc - 3);
    get_vector(&x, argc - 2, -1);
    get_vector(&y, argc - 3, -1);
    drop = FALSE;
  } else if (argc == 4) {
    y_ref = yget_ref(argc - 3); /* get ref *before* yarg_typeid */
    if (yarg_typeid(argc - 3) == Y_RANGE) {
      y_ref = yget_ref(argc - 4); /* get ref *before* get_vector */
      get_vector(&x, argc - 2, argc - 3);
      get_vector(&y, argc - 4, -1);
      drop = FALSE;
    } else {
      get_vector(&x, argc - 2, -1);
      get_vector(&y, argc - 3, argc - 4);
      drop = TRUE;
    }
  } else if (argc == 5) {
    y_ref = yget_ref(argc - 4); /* get ref *before* get_vector */
    get_vector(&x, argc - 2, argc - 3);
    get_vector(&y, argc - 4, argc - 5);
    drop = TRUE;
  } else {
    y_error("expecting 3, 4 or 5 arguments");
    return;
  }
  if (y_ref < 0) {
    y_error("argument Y must be a variable (not an expression)");
    return;
  }

  /* Check sizes and types. */
  n = check_vector_sizes(&x, &y);
  type = Y_FLOAT;
  if (type < VECTOR_TYPE(x)) {
    type = VECTOR_TYPE(x);
  }
  if (type < VECTOR_TYPE(y)) {
    type = VECTOR_TYPE(y);
  }
  if (VECTOR_TYPE(x) != type) {
    coerce_vector(&x, type);
  }
  if (VECTOR_TYPE(y) != type) {
    coerce_vector(&y, type);
  } else {
    y_ref = -1L; /* no needs to re-define global symbol */
  }
  if (type == Y_DOUBLE) {
    double alpha = ygets_d(argc - 1);
    CALL_DAXPY(n,alpha,x, y);
  } else if (type == Y_FLOAT) {
    float alpha = ygets_f(argc - 1);
    CALL_SAXPY(n,alpha,x,y);
  } else {
    double alpha[2];
    get_z(argc - 1, alpha);
    CALL_ZAXPY(n,alpha,x,y);
  }

  /* Make sure Y is updated and that it is left on top of the stack. */
  if (y_ref >= 0L) {
    yput_global(y_ref, y.arr.iarg);
  }
  if (drop) {
    yarg_drop(1);
  }
}

/*---------------------------------------------------------------------------*/
/* LEVEL 2 BLAS */

void Y_lpk_gemv(int argc)
{
  matrix_t a;
  vector_t x, y;
  long y_ref;
#ifdef USE_CBLAS
  CBLAS_TRANSPOSE trans;
#else
  CHARACTER trans[2];
#endif
  BLAS_INTEGER m, n;
  int type, i;
  int alpha_iarg, beta_iarg;

  /* Get arguments (not ALPHA nor BETA). */
  if (argc < 6 || argc > 8) {
    y_error("lpk_gemv takes 6 to 8 arguments");
  }
  GET_BLAS_TRANS(argc - 1, trans);
  alpha_iarg = argc - 2;
  get_matrix(&a, argc - 3);
  if (argc == 8 || (argc == 7 && yarg_typeid(argc - 5) == Y_RANGE)) {
    get_vector(&x, argc - 4, argc - 5);
    beta_iarg = argc - 6;
  } else {
    get_vector(&x, argc - 4, -1);
    beta_iarg = argc - 5;
  }
  y_ref = yget_ref(beta_iarg - 1);
  if (y_ref < 0L) {
    y_error("argument Y must be a variable");
  }
  get_vector(&y, beta_iarg - 1, beta_iarg - 2);

  /* Check dimensions. */
  if (MATRIX_RANK(a) != x.rank + y.rank) {
    goto badDims;
  }
  if (BLAS_ALTERNATIVE(trans != CblasNoTrans, trans[0] != 'N')) {
    /* Check that A is XDIMS-by-YDIMS. */
    if (x.rank < 2) {
      if (a.arr.dims[1] != x.size) {
        goto badDims;
      }
    } else {
      for (i = 1; i <= x.rank; ++i) {
        if (a.arr.dims[i] != x.arr.dims[i]) {
          goto badDims;
        }
      }
    }
    if (y.rank < 2) {
      if (a.arr.dims[MATRIX_RANK(a)] != y.size) {
        goto badDims;
      }
    } else {
      for (i = 1; i <= y.rank; ++i) {
        if (a.arr.dims[x.rank + i] != y.arr.dims[i]) {
          goto badDims;
        }
      }
    }
    a.nrows = x.size;
    a.ncols = y.size;
  } else {
    /* Check that A is YDIMS-by-XDIMS. */
    if (y.rank < 2) {
      if (a.arr.dims[1] != y.size) {
        goto badDims;
      }
    } else {
      for (i = 1; i <= y.rank; ++i) {
        if (a.arr.dims[i] != y.arr.dims[i]) {
          goto badDims;
        }
      }
    }
    if (x.rank < 2) {
      if (a.arr.dims[MATRIX_RANK(a)] != x.size) {
        goto badDims;
      }
    } else {
      for (i = 1; i <= x.rank; ++i) {
        if (a.arr.dims[y.rank + i] != x.arr.dims[i]) {
          goto badDims;
        }
      }
    }
    a.nrows = y.size;
    a.ncols = x.size;
  }

  /* Figure out the data type and convert arguments. */
  type = Y_FLOAT;
  if (type < MATRIX_TYPE(a)) {
    type = MATRIX_TYPE(a);
  }
  if (type < VECTOR_TYPE(x)) {
    type = VECTOR_TYPE(x);
  }
  if (type < VECTOR_TYPE(y)) {
    type = VECTOR_TYPE(y);
  }
  if (VECTOR_TYPE(x) != type) {
    coerce_vector(&x, type);
  }
  if (VECTOR_TYPE(y) != type) {
    coerce_vector(&y, type);
  } else {
    y_ref = -1L; /* no needs to re-define global symbol */
  }
  if (MATRIX_TYPE(a) != type) {
    coerce_matrix(&a, type);
  }

  /* Apply the operation. */
  m = a.nrows;
  n = a.ncols;
  if (type == Y_DOUBLE) {
    double alpha = ygets_d(alpha_iarg);
    double beta = ygets_d(beta_iarg);
    CALL_DGEMV(trans, m, n, alpha, a, x, beta, y);
  } else if (type == Y_FLOAT) {
    float alpha = ygets_f(alpha_iarg);
    float beta = ygets_f(beta_iarg);
    CALL_SGEMV(trans, m, n, alpha, a, x, beta, y);
  } else {
    double alpha[2], beta[2];
    get_z(alpha_iarg, alpha);
    get_z(beta_iarg, beta);
    CALL_ZGEMV(trans, m, n, alpha, a, x, beta, y);
  }

  if (y_ref >= 0L) {
    yput_global(y_ref, y.arr.iarg);
  }
  if (beta_iarg == 2) {
    yarg_drop(1);
  }
  return;
 badDims:
  y_error("incompatible dimensions");
}

/*---------------------------------------------------------------------------*/
/* LEVEL 3 BLAS */

/* C := alpha*op( A )*op( B ) + beta*C */
void Y_lpk_gemm(int argc)
{
  matrix_t a, b, c;
  long c_ref; /* index to variable C */
  long k_size, m_size, n_size; /* number of elements */
  long k_rank, m_rank, n_rank; /* rank */
  long i, a_off, b_off; /* index and offsets along dimension lists */
  long len; /* length of dimension */
#ifdef USE_CBLAS
  CBLAS_TRANSPOSE transa, transb;
#else
  CHARACTER transa[2], transb[2];
#endif
  BLAS_INTEGER k, m, n;
  int type;
  int a_rev, b_rev; /* reverse dimensions of A or B? */

  /* Parse arguments (ALPHA and BETA are extracted later). */
  if (argc != 7) {
    y_error("blas_gemm takes 7 arguments");
  }
  GET_BLAS_TRANS(argc - 1, transa);
  GET_BLAS_TRANS(argc - 2, transb);
  get_matrix(&a, argc - 4);
  get_matrix(&b, argc - 5);
  c_ref = yget_ref(argc - 7);
  if (c_ref < 0L) {
    y_error("argument C must be a variable");
  }
  get_matrix(&c, argc - 7);

  /* Check number of dimensions. */
  k_rank = MATRIX_RANK(a) + MATRIX_RANK(b) - MATRIX_RANK(c);
  m_rank = MATRIX_RANK(a) + MATRIX_RANK(c) - MATRIX_RANK(b);
  n_rank = MATRIX_RANK(b) + MATRIX_RANK(c) - MATRIX_RANK(a);
  if (k_rank < 0L || m_rank < 0L || n_rank < 0L || (k_rank & 1L) != 0L) {
  incompatible:
    y_error("arguments A, B and C have incompatible dimensions");
  }
  k_rank >>= 1;
  m_rank >>= 1;
  n_rank >>= 1;

  /* Check dimensions and count elements along dimensions M, K and N. */
  a_rev = BLAS_ALTERNATIVE(transa != CblasNoTrans, transa[0] != 'N');
  b_rev = BLAS_ALTERNATIVE(transb != CblasNoTrans, transb[0] != 'N');
  m_size = 1L;
  a_off = (a_rev ? k_rank : 0L);
  for (i = 1L; i <= m_rank; ++i) {
    if ((len = a.arr.dims[a_off + i]) != c.arr.dims[i])  {
      goto incompatible;
    }
    m_size *= len;
  }
  n_size = 1L;
  b_off = (b_rev ? 0L : k_rank);
  for (i = 1L; i <= n_rank; ++i) {
    if ((len = b.arr.dims[b_off + i]) != c.arr.dims[m_rank + i])  {
      goto incompatible;
    }
    n_size *= len;
  }
  k_size = 1L;
  a_off = (a_rev ? 0L : m_rank);
  b_off = (b_rev ? n_rank : 0L);
  for (i = 1L; i <= k_rank; ++i) {
    if ((len = a.arr.dims[a_off + i]) != b.arr.dims[b_off + i])  {
      goto incompatible;
    }
    k_size *= len;
  }
  if (a_rev) {
    a.nrows = k_size;
    a.ncols = m_size;
  } else {
    a.nrows = m_size;
    a.ncols = k_size;
  }
  if (b_rev) {
    b.nrows = n_size;
    b.ncols = k_size;
  } else {
    b.nrows = k_size;
    b.ncols = n_size;
  }
  c.nrows = m_size;
  c.ncols = n_size;

  /* Figure out the data type and convert arguments. */
  type = Y_FLOAT;
  if (type < MATRIX_TYPE(a)) {
    type = MATRIX_TYPE(a);
  }
  if (type < MATRIX_TYPE(b)) {
    type = MATRIX_TYPE(b);
  }
  if (type < MATRIX_TYPE(c)) {
    type = MATRIX_TYPE(c);
  }
  if (MATRIX_TYPE(a) != type) {
    coerce_matrix(&a, type);
  }
  if (MATRIX_TYPE(b) != type) {
    coerce_matrix(&b, type);
  }
  if (MATRIX_TYPE(c) != type) {
    coerce_matrix(&c, type);
  } else {
    c_ref = -1L; /* no needs to re-define global symbol */
  }

  /* Apply the operation. */
  k = k_size;
  m = m_size;
  n = n_size;
  if (type == Y_DOUBLE) {
    double alpha = ygets_d(argc - 3);
    double beta = ygets_d(argc - 6);
    CALL_DGEMM(transa, transb, m, n, k, alpha, a, b, beta, c);
  } else if (type == Y_FLOAT) {
    float alpha = ygets_f(argc - 3);
    float beta = ygets_f(argc - 6);
    CALL_SGEMM(transa, transb, m, n, k, alpha, a, b, beta, c);
  } else {
    double alpha[2], beta[2];
    get_z(argc - 3, alpha);
    get_z(argc - 6, beta);
    CALL_ZGEMM(transa, transb, m, n, k, alpha, a, b, beta, c);
  }
  if (c_ref >= 0L) {
    yput_global(c_ref, c.arr.iarg);
  }
}

void Y_lpk_syrk(int argc)
{
  matrix_t a, c;
  long c_ref, dim;
  BLAS_INTEGER n, k;
#ifdef USE_CBLAS
  CBLAS_UPLO uplo;
  CBLAS_TRANSPOSE trans;
#else
  CHARACTER uplo[2], trans[2];
#endif
  int c_copy, type, i, p, r;

  if (argc != 6) {
    y_error("blas_syrk takes 6 arguments");
  }
  c_ref = yget_ref(argc - 6);
  if (yarg_subroutine()) {
    /* operation is done in-place */
    c_copy = FALSE;
    if (c_ref < 0) {
      y_error("argument C must be a variable in this context");
      return;
    }
  } else {
    /* no needs to copy if C is an expression */
    c_copy = (c_ref >= 0L);
  }
  GET_BLAS_UPLO(argc - 1, uplo);
  GET_BLAS_TRANS(argc - 2, trans);
  get_matrix(&a, argc - 4);
  get_matrix(&c, argc - 6);

  /* Check dimensions and data type. */
  if (MATRIX_RANK(c) % 2 != 0) {
      goto badDimsForC;
  }
  p = MATRIX_RANK(c)/2;
  c.nrows = 1L;
  for (i = 1; i <= p; ++i) {
    if ((dim = c.arr.dims[i]) != c.arr.dims[p + i]) {
      goto badDimsForC;
    }
    c.nrows *= dim;
  }
  c.ncols = c.nrows;
  n = MATRIX_NROWS(c);
  if (MATRIX_RANK(a) < p) {
    goto incompatibleDims;
  }
  if (BLAS_ALTERNATIVE(trans == CblasNoTrans, trans[0] == 'N')) {
    /* A must be N-by-K. */
    for (i = 1; i <= p; ++i) {
      if (a.arr.dims[i] != c.arr.dims[i]) {
        goto incompatibleDims;
      }
    }
    a.nrows = c.nrows;
    a.ncols = 1L;
    for (i = p + 1; i <= MATRIX_RANK(a); ++i) {
      a.ncols *= a.arr.dims[i];
    }
    k = MATRIX_NCOLS(a);
  } else {
    /* A must be K-by-N. */
    r = MATRIX_RANK(a) - p;
    for (i = 1; i <= p; ++i) {
      if (a.arr.dims[r + i] != c.arr.dims[i]) {
        goto incompatibleDims;
      }
    }
    a.ncols = c.nrows;
    a.nrows = 1L;
    for (i = 1; i <= r; ++i) {
      a.nrows *= a.arr.dims[i];
    }
    k = MATRIX_NROWS(a);
  }
  type = Y_FLOAT;
  if (type < MATRIX_TYPE(a)) {
    type = MATRIX_TYPE(a);
  }
  if (type < MATRIX_TYPE(c)) {
    type = MATRIX_TYPE(c);
  }
  if (type == Y_COMPLEX) {
    y_error("arguments must not be of type complex");
    return;
  }

  /* Convert arguments. */
  if (MATRIX_TYPE(a) != type) {
    coerce_matrix(&a, type);
  }
  if (MATRIX_TYPE(c) != type) {
    coerce_matrix(&c, type);
    if (c_ref >= 0L) {
      yput_global(c_ref, c.arr.iarg);
    }
  } else if (c_copy) {
    copy_array(&c.arr);
  }

  /* Apply the operation. */
  if (type == Y_DOUBLE) {
    double alpha = ygets_d(argc - 3);
    double beta = ygets_d(argc - 5);
    CALL_DSYRK(uplo, trans, n, k, alpha, a, beta, c);
  } else {
    float alpha = ygets_f(argc - 3);
    float beta = ygets_f(argc - 5);
    CALL_SSYRK(uplo, trans, n, k, alpha, a, beta, c);
  }
  return;
 incompatibleDims:
  y_error("matrices A and C have incompatible dimensions");
  return;
 badDimsForC:
  y_error("argument C must be a square matrix");
  return;
}


/*---------------------------------------------------------------------------*/
/* LAPACK */

void Y_lpk_potrf(int argc)
{
  long ntot, ref;
  long dims[Y_DIMSIZE];
  void *a;
  CHARACTER uplo[2];
  INTEGER n, info;
  int copy, type, save, i, r;

  if (argc != 2) {
    y_error("lpk_potrf takes 2 arguments");
  }
  get_lapack_uplo(argc - 1, uplo);
  ref = yget_ref(argc - 2);
  if (yarg_subroutine()) {
    if (ref < 0) {
      y_error("argument A must be a variable in this context");
    }
    copy = FALSE;
    save = TRUE;
  } else {
    copy = (ref >= 0L);
    save = FALSE;
  }
  type = yarg_typeid(argc - 2);
  if (type == Y_DOUBLE) {
    a = ygeta_d(argc - 2, &ntot, dims);
    save = FALSE;
  } else if (type <= Y_FLOAT) {
    a = ygeta_f(argc - 2, &ntot, dims);
    if (type != Y_FLOAT) {
      type = Y_FLOAT;
      copy = FALSE;
    } else {
      save = FALSE;
    }
  } else {
    y_error("unsupported matrix type");
    return;
  }
  r = dims[0];
  if (r % 2 != 0) {
    goto badDimList;
  }
  r >>= 1;
  n = 1;
  for (i = 1; i <= r; ++i) {
    if (dims[i] != dims[r + i]) {
      goto badDimList;
    }
    n *= dims[i];
  }
  if (type == Y_DOUBLE) {
    if (copy) {
      /* FIXME: just copy the relevant triangular part */
      a = push_copy_d(a, dims);
    }
    DPOTRF(uplo, &n, a, &n, &info);
  } else {
    if (copy) {
      /* FIXME: just copy the relevant triangular part */
      a = push_copy_f(a, dims);
    }
    SPOTRF(uplo, &n, a, &n, &info);
  }
  if (info != 0) {
    char msg[80];
    if (info < 0) {
      sprintf(msg, "illegal value for %d-th argument in POTRF", -info);
    } else {
      sprintf(msg, "the leading minor of order %d is not positive definite",
             info);
    }
    y_error(msg);
    return;
  }
  if (save && ref >= 0L) {
    yput_global(ref, 0);
  }
  return;

 badDimList:
  y_error("argument A must be a N-by-N matrix");
  return;
}

void Y_lpk_potrs(int argc)
{
  linear_system_t sys;
  CHARACTER uplo[2];
  INTEGER info;

  if (argc != 3) {
    y_error("lpk_potrs takes three arguments");
  }

  /* uplo flag */
  get_lapack_uplo(2, uplo);

  /* get LHS matrix and RHS vector */
  get_linear_system(&sys, 1, 0, REAL_SYSTEM);

  /* solve the linear system */
  if (sys.type == Y_DOUBLE) {
    (void)DPOTRS(uplo, &sys.n, &sys.nrhs, sys.a, &sys.lda,
                 sys.b, &sys.ldb, &info);
  } else /* Y_FLOAT */ {
    (void)SPOTRS(uplo, &sys.n, &sys.nrhs, sys.a, &sys.lda,
                 sys.b, &sys.ldb, &info);
  }
  if (info != 0) {
    y_error(info > 0 ? "singular matrix" : "bug in lpk_potrs");
  }
}

/* FIXME: merge this with Y_lpk_potrf */
/* FIXME: copy the other triangular part of the matrix? */
/* FIXME: if not in-place, just copy the relevant triangular part */
/* FIXME: same remark for data conversion */
void Y_lpk_potri(int argc)
{
  long ntot, ref;
  long dims[Y_DIMSIZE];
  void *a;
  CHARACTER uplo[2];
  INTEGER n, info;
  int copy, type, save, i, r;

  if (argc != 2) {
    y_error("lpk_potri takes two arguments");
  }

  /* uplo flag */
  get_lapack_uplo(argc - 1, uplo);
  ref = yget_ref(argc - 2);
  if (yarg_subroutine()) {
    if (ref < 0) {
      y_error("argument A must be a variable in this context");
    }
    copy = FALSE;
    save = TRUE;
  } else {
    copy = (ref >= 0L);
    save = FALSE;
  }
  type = yarg_typeid(argc - 2);
  if (type == Y_DOUBLE) {
    a = ygeta_d(argc - 2, &ntot, dims);
    save = FALSE;
  } else if (type <= Y_FLOAT) {
    a = ygeta_f(argc - 2, &ntot, dims);
    if (type != Y_FLOAT) {
      type = Y_FLOAT;
      copy = FALSE;
    } else {
      save = FALSE;
    }
  } else {
    y_error("unsupported matrix type");
    return;
  }
  r = dims[0];
  if (r % 2 != 0) {
    goto badDimList;
  }
  r >>= 1;
  n = 1;
  for (i = 1; i <= r; ++i) {
    if (dims[i] != dims[r + i]) {
      goto badDimList;
    }
    n *= dims[i];
  }
  if (type == Y_DOUBLE) {
    if (copy) {
      /* FIXME: just copy the relevant triangular part */
      a = push_copy_d(a, dims);
    }
    DPOTRI(uplo, &n, a, &n, &info);
  } else {
    if (copy) {
      /* FIXME: just copy the relevant triangular part */
      a = push_copy_f(a, dims);
    }
    SPOTRI(uplo, &n, a, &n, &info);
  }
  if (info != 0) {
    char msg[256];
    if (info < 0) {
      sprintf(msg, "illegal value for %d-th argument in POTRI", -info);
    } else {
      sprintf(msg, "the (%d,%d) element of the factor U or L is zero, " \
              "and the inverse could not be computed", info, info);
    }
    y_error(msg);
    return;
  }
  if (save && ref >= 0L) {
    yput_global(ref, 0);
  }
  return;

 badDimList:
  y_error("argument A must be a N-by-N matrix");
  return;
}

void Y_lpk_gesv(int argc)
{
  linear_system_t sys;
  INTEGER *ipvt;
  INTEGER info;
  long n;

  if (argc != 2) {
    y_error("lpk_gesv takes two arguments");
  }

  /* get LHS matrix and RHS vector */
  n = get_linear_system(&sys, 1, 0, REAL_OR_COMPLEX_SYSTEM);

  /* push temporary workspace on top of the stack */
  ipvt = (INTEGER *)ypush_scratch(sizeof(INTEGER)*n, NULL);

  /* solve the linear system */
  if (sys.type == Y_DOUBLE) {
    (void)DGESV(&sys.n, &sys.nrhs, sys.a, &sys.lda, ipvt,
                sys.b, &sys.ldb, &info);
  } else if (sys.type == Y_FLOAT) {
    (void)SGESV(&sys.n, &sys.nrhs, sys.a, &sys.lda, ipvt,
                sys.b, &sys.ldb, &info);
  } else /* Y_COMPLEX */ {
    (void)ZGESV(&sys.n, &sys.nrhs, sys.a, &sys.lda, ipvt,
                sys.b, &sys.ldb, &info);
  }
  if (info != 0) {
    y_error(info > 0 ? "singular matrix" : "bug in lpk_gesv");
  }

  /* get rid of temporary workspace on top of the stack */
  yarg_drop(1);
}

void Y_lpk_posv(int argc)
{
  linear_system_t sys;
  CHARACTER uplo[2];
  INTEGER info;

  if (argc != 3) {
    y_error("lpk_posv takes three arguments");
  }

  /* uplo flag */
  get_lapack_uplo(2, uplo);

  /* get LHS matrix and RHS vector */
  get_linear_system(&sys, 1, 0, REAL_OR_COMPLEX_SYSTEM);

  /* solve the linear system */
  if (sys.type == Y_DOUBLE) {
    (void)DPOSV(uplo, &sys.n, &sys.nrhs, sys.a, &sys.lda,
                sys.b, &sys.ldb, &info);
  } else if (sys.type == Y_FLOAT) {
    (void)SPOSV(uplo, &sys.n, &sys.nrhs, sys.a, &sys.lda,
                sys.b, &sys.ldb, &info);
  } else /* Y_COMPLEX */ {
    (void)ZPOSV(uplo, &sys.n, &sys.nrhs, sys.a, &sys.lda,
                sys.b, &sys.ldb, &info);
  }
  if (info != 0) {
    y_error(info > 0 ? "singular matrix" : "bug in lpk_posv");
  }
}

static void sysv_or_hesv(int argc, int which)
{

  linear_system_t sys;
  CHARACTER uplo[2];
  INTEGER info, *ipvt, lwork;
  void *work;
  long n;

  if (argc != 3) {
    y_error((which == 0 ? "lpk_sysv takes three arguments" :
             "lpk_hesv takes three arguments"));
  }

  /* uplo flag */
  get_lapack_uplo(2, uplo);

  /* get LHS matrix and RHS vector */
  n = get_linear_system(&sys, 1, 0,
                        (which == 0 ? REAL_OR_COMPLEX_SYSTEM
                         : COMPLEX_SYSTEM));

  /* push temporary workspace on top of the stack */
  ipvt = (INTEGER *)ypush_scratch(sizeof(INTEGER)*n, NULL);

  /* solve the linear system */
  if (sys.type == Y_DOUBLE) {
    double dummy;
    work = &dummy;
    lwork = -1;
    (void)DSYSV(uplo, &sys.n, &sys.nrhs, sys.a, &sys.lda, ipvt,
                sys.b, &sys.ldb, work, &lwork, &info);
    if (info == 0) {
      if ((lwork = (INTEGER)dummy) > 1) {
        work = ypush_scratch(sizeof(double)*lwork, NULL);
      } else {
        lwork = 1;
      }
      (void)DSYSV(uplo, &sys.n, &sys.nrhs, sys.a, &sys.lda, ipvt,
                  sys.b, &sys.ldb, work, &lwork, &info);
    }
  } else if (sys.type == Y_FLOAT) {
    float dummy;
    work = &dummy;
    lwork = -1;
    (void)SSYSV(uplo, &sys.n, &sys.nrhs, sys.a, &sys.lda, ipvt,
                sys.b, &sys.ldb, work, &lwork, &info);
    if (info == 0) {
      if ((lwork = (INTEGER)dummy) > 1) {
        work = ypush_scratch(sizeof(float)*lwork, NULL);
      } else {
        lwork = 1;
      }
      (void)SSYSV(uplo, &sys.n, &sys.nrhs, sys.a, &sys.lda, ipvt,
                  sys.b, &sys.ldb, work, &lwork, &info);
    }
  } else /* Y_COMPLEX */ {
    double dummy[2];
    work = dummy;
    lwork = -1;
    if (which == 0) {
      (void)ZSYSV(uplo, &sys.n, &sys.nrhs, sys.a, &sys.lda, ipvt,
                  sys.b, &sys.ldb, work, &lwork, &info);
    } else {
      (void)ZHESV(uplo, &sys.n, &sys.nrhs, sys.a, &sys.lda, ipvt,
                  sys.b, &sys.ldb, work, &lwork, &info);
    }
    if (info == 0) {
      if ((lwork = (INTEGER)dummy[0]) > 1) {
        work = ypush_scratch((2*sizeof(double))*lwork, NULL);
      } else {
        lwork = 1;
      }
      if (which == 0) {
        (void)ZSYSV(uplo, &sys.n, &sys.nrhs, sys.a, &sys.lda, ipvt,
                    sys.b, &sys.ldb, work, &lwork, &info);
      } else {
        (void)ZHESV(uplo, &sys.n, &sys.nrhs, sys.a, &sys.lda, ipvt,
                    sys.b, &sys.ldb, work, &lwork, &info);
      }
    }
  }
  if (info != 0) {
    y_error(info > 0 ? "singular matrix" :
            (which == 0 ? "bug in lpk_sysv" : "bug in lpk_hesv"));
  }

  /* get rid of temporary workspaces on top of the stack */
  yarg_drop((lwork > 1 ? 2 : 1));
}

void Y_lpk_sysv(int argc)
{
  sysv_or_hesv(argc, 0);
}

void Y_lpk_hesv(int argc)
{
  sysv_or_hesv(argc, 1);
}

#if 0 /* not used */
static void free_global(long index)
{
  if (index >= 0L) {
    ypush_nil();
    yput_global(index, 0);
    yarg_drop(1);
  }
}
#endif

/* Force a numerical array stored on the stack to be a temporary one and
   returns the address of the array. */
static void *force_scratch(int iarg, void *ptr, int type,
                           long ntot, long *dims);

static void *push_1d(int type, long n);
static void *push_2d(int type, long m, long n);

/*---------------------------------------------------------------------------*/
/* EIGENVALUES AND EIGENVECTORS */

void Y_lpk_eigen(int argc)
{
  CHARACTER uplo[2], jobz[2];
  long dims[Y_DIMSIZE];
  void *a, *w, *work, *rwork, *iwork;
  long index, v_ref, ntot, offset;
  INTEGER n, info, lda, lwork, lrwork, liwork;
  int slow, iarg, type, positional_args;

  /* Parse arguments. */
  slow = FALSE;
  a = NULL;
  v_ref = -1L;
  type = Y_VOID;
  n = 0;
  positional_args = 0;
  for (iarg = argc - 1; iarg >= 0; --iarg) {
    index = yarg_key(iarg);
    if (index >= 0L) {
      --iarg;
      if (index == slow_index) {
        slow = yarg_true(iarg);
      } else {
        y_error("unknown keyword");
      }
    } else {
      ++positional_args;
      if (positional_args == 1) {
        /* Get UPLO flag. */
        get_lapack_uplo(iarg, uplo);
      } else if (positional_args == 2) {
        /* Get argument A making sure it is a temporary array because SYEV(D)
           and HEEV(D) destroy the input matrix. */
        a = ygeta_any(iarg, &ntot, dims, &type);
        if (type > Y_COMPLEX || dims[0] !=  2 || dims[1] != dims[2]) {
          y_error("expecting N-by-N real or complex matrix A");
        }
        if (type < Y_FLOAT) {
          a = ygeta_coerce(iarg, a, ntot, dims, type, Y_FLOAT);
          type = Y_FLOAT;
        } else {
          a = force_scratch(iarg, a, type, ntot, dims);
        }
        n = dims[1];
      } else if (positional_args == 3) {
        /* Get argument V. */
        v_ref = yget_ref(iarg);
      } else {
        y_error("too many arguments");
      }
    }
  }
  if (positional_args < 1) {
    y_error("too few arguments");
  }

  /* The routine will push 2 items on top of the stack: the array W to store
     eigenvalues and a workspace. (This is less than 8 items, so there are no
     needs to call ypush_check in principle.) */
  ypush_check(2);
  
  /* Create W on top of the stack. */
  w = push_1d((type == Y_FLOAT ? Y_FLOAT : Y_DOUBLE), n);

  /* Perform the decomposition.  First get the optimal LWORK size, then
     allocate workspaces and execute the operations. */
  lda = n;
  lwork = -1;
  liwork = -1;
  info = 0;
  jobz[0] = (v_ref < 0L ? 'N' : 'V');
  jobz[1] = '\0';
  if (slow) {
    /* Simple version of the algorithm.*/
    if (type == Y_FLOAT) {
      float temp;
      SSYEV(jobz, uplo, &n, a, &lda, w, &temp, &lwork, &info);
      if (info == 0) {
        lwork = (INTEGER)temp;
        work = ypush_scratch(lwork*sizeof(float), NULL);
        SSYEV(jobz, uplo, &n, a, &lda, w, work, &lwork, &info);
      }
    } else if (type == Y_DOUBLE) {
      double temp;
      DSYEV(jobz, uplo, &n, a, &lda, w, &temp, &lwork, &info);
      if (info == 0) {
        lwork = (INTEGER)temp;
        work = ypush_scratch(lwork*sizeof(double), NULL);
        DSYEV(jobz, uplo, &n, a, &lda, w, work, &lwork, &info);
      }
    } else {
      double temp[2];
      ZHEEV(jobz, uplo, &n, a, &lda, w, temp, &lwork, NULL, &info);
      if (info == 0) {
        lrwork = 3*MAX(n,1) - 2;
        lwork = (INTEGER)temp[0];
        work = ypush_scratch((2*lwork + lrwork)*sizeof(double), NULL);
        rwork = INCR_ADDRESS(work, 2*lwork*sizeof(double));
        ZHEEV(jobz, uplo, &n, a, &lda, w, work, &lwork, rwork, &info);
      }
    }
  } else {
    /* Divide-and-conquer version of the algorithm. */
    if (type == Y_FLOAT) {
      float temp;
      INTEGER itemp;
      SSYEVD(jobz, uplo, &n, a, &lda, w, &temp, &lwork,
             &itemp, &liwork, &info);
      if (info == 0) {
        lwork = (INTEGER)temp;
        liwork = itemp;
        offset = ROUND_UP(lwork*sizeof(float), sizeof(INTEGER));
        work = ypush_scratch(offset + liwork*sizeof(INTEGER), NULL);
        iwork = INCR_ADDRESS(work, offset);
        SSYEVD(jobz, uplo, &n, a, &lda, w, work, &lwork,
               iwork, &liwork, &info);
      }
    } else if (type == Y_DOUBLE) {
      double temp;
      INTEGER itemp;
      DSYEVD(jobz, uplo, &n, a, &lda, w, &temp, &lwork,
             &itemp, &liwork, &info);
      if (info == 0) {
        lwork = (INTEGER)temp;
        liwork = itemp;
        offset = ROUND_UP(lwork*sizeof(double), sizeof(INTEGER));
        work = ypush_scratch(offset + liwork*sizeof(INTEGER), NULL);
        iwork = INCR_ADDRESS(work, offset);
        DSYEVD(jobz, uplo, &n, a, &lda, w, work, &lwork,
               iwork, &liwork, &info);
      }
    } else {
      double temp[2], rtemp;
      INTEGER itemp;
      ZHEEVD(jobz, uplo, &n, a, &lda, w, temp, &lwork, &rtemp, &lrwork,
             &itemp, &liwork, &info);
      if (info == 0) {
        lwork = (INTEGER)temp[0];
        liwork = itemp;
        lrwork = (INTEGER)rtemp;
        offset = ROUND_UP((lrwork + 2*lwork)*sizeof(double), sizeof(INTEGER));
        work = ypush_scratch(offset + liwork*sizeof(INTEGER), NULL);
        rwork = INCR_ADDRESS(work, 2*lwork*sizeof(double));
        iwork = INCR_ADDRESS(work, offset);
        ZHEEVD(jobz, uplo, &n, a, &lda, w, work, &lwork, rwork, &lrwork,
               iwork, &liwork, &info);
      }
    }
  }

  /* Build result. */
  if (info != 0) {
    /* Algorithm failed.  raise an error if called as a subroutine. */
    char *msg;
    if (info > 0) {
      msg = "algorithm did not converge";
    } else {
        sprintf(msgbuf,
                "bug in x%s%s wrapper: %d-th argument had an illegal value",
                (type == Y_COMPLEX ? "HEEV" : "SYEV" ),
                (slow ? "" : "D"), -(int)info);
        msg = msgbuf;
    }
    y_error(msg);
  }
  yarg_drop(1); /* drop workspace */
  if (v_ref >= 0L) {
    yput_global(v_ref, 1);
  }
}

/*---------------------------------------------------------------------------*/
/* SINGULAR VALUE DECOMPOSITION */

static void gesvd(int argc, int divide_and_conquer);

void Y_lpk_gesvd(int argc)
{
  gesvd(argc, 0);
}

void Y_lpk_gesdd(int argc)
{
  gesvd(argc, 1);
}

static void gesvd(int argc, int divide_and_conquer)
{
  long dims[Y_DIMSIZE];
  void *a, *s, *u, *vt, *work, *rwork, *iwork;
  long index, s_ref, u_ref, vt_ref, ntot, offset, size;
  INTEGER m, n, info, lda, ldu, ldvt, lwork, lrwork;
  int full, iarg, type, positional_args;

  /* Parse arguments. */
  full = FALSE;
  a = NULL;
  s_ref = -1L;
  u_ref = -1L;
  vt_ref = -1L;
  type = Y_VOID;
  m = 0;
  n = 0;
  positional_args = 0;
  for (iarg = argc - 1; iarg >= 0; --iarg) {
    index = yarg_key(iarg);
    if (index >= 0L) {
      --iarg;
      if (index == full_index) {
        full = yarg_true(iarg);
      } else {
        y_error("unknown keyword");
      }
    } else {
      ++positional_args;
      if (positional_args == 1) {
        /* Get argument A making sure it is a temporary array because GESVD
           and GESDD destroy the input matrix. */
        type = yarg_typeid(iarg);
        if (type > Y_COMPLEX || yarg_rank(iarg) !=  2) {
          y_error("expecting a real or complex matrix");
        }
        if (type <= Y_FLOAT) {
          a = ygeta_f(iarg, &ntot, dims);
          type = Y_FLOAT;
        } else if (type == Y_DOUBLE) {
          a = ygeta_d(iarg, &ntot, dims);
        } else {
          a = ygeta_z(iarg, &ntot, dims);
        }
        a = force_scratch(iarg, a, type, ntot, dims);
        m = dims[1];
        n = dims[2];
      } else if (positional_args == 2) {
        /* Get argument S. */
        s_ref = yget_ref(iarg);
      } else if (positional_args == 3) {
        /* Get argument U. */
        u_ref = yget_ref(iarg);
      } else if (positional_args == 4) {
        /* Get argument VT. */
        vt_ref = yget_ref(iarg);
      } else {
        y_error("too many arguments");
      }
    }
  }
  if (positional_args < 1) {
    y_error("too few arguments");
  }

  /* The routine will push 4 items on top of the stack: S, U, VT and a
     workspace. (This is less than 8 items, so there are no needs to call
     ypush_check in principle.) */
  ypush_check(4);
  
  /* Create S, U and VT on top of the stack.  FIXME: use A to save storage. */
  s = push_1d((type == Y_FLOAT ? Y_FLOAT : Y_DOUBLE), MIN(m, n));
  if (u_ref >= 0L || (divide_and_conquer && vt_ref >= 0L)) {
    ldu = m;
    u = push_2d(type, ldu, (full ? m : MIN(m, n)));
  } else {
    ldu = 1;
    u = NULL;
    ypush_nil();
  }
  if (vt_ref >= 0L || (divide_and_conquer && u_ref >= 0L)) {
    ldvt = (full ? n : MIN(m, n));
    vt = push_2d(type, ldvt, n);
  } else {
    ldvt = 1;
    vt = NULL;
    ypush_nil();
  }

  /* Perform the decomposition.  First get the optimal LWORK size, then
     allocate workspaces and execute the operations. */
  lda = m;
  lwork = -1;
  info = 0;
  if (divide_and_conquer) {
    /* Divide-and-conquer version of the algorithm. */
    CHARACTER jobz[2];
    jobz[0] = (u == NULL ? 'N' : (full ? 'A' : 'S'));
    jobz[1] = '\0';
    if (type == Y_FLOAT) {
      float temp;
      SGESDD(jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
             &temp, &lwork, NULL, &info);
      if (info == 0) {
        lwork = (INTEGER)temp;
        offset = ROUND_UP(lwork*sizeof(float), sizeof(INTEGER));
        size = offset + 8*MIN(m,n)*sizeof(INTEGER);
        work = ypush_scratch(size, NULL);
        iwork = INCR_ADDRESS(work, offset);
        SGESDD(jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
               work, &lwork, iwork, &info);
      }
    } else if (type == Y_DOUBLE) {
      double temp;
      DGESDD(jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
             &temp, &lwork, NULL, &info);
      if (info == 0) {
        lwork = (INTEGER)temp;
        offset = ROUND_UP(lwork*sizeof(double), sizeof(INTEGER));
        size = offset + 8*MIN(m,n)*sizeof(INTEGER);
        work = ypush_scratch(size, NULL);
        iwork = INCR_ADDRESS(work, offset);
        DGESDD(jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
               work, &lwork, iwork, &info);
      }
    } else {
      double temp[2];
      ZGESDD(jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
             temp, &lwork, NULL, NULL, &info);
      if (info == 0) {
        lwork = (INTEGER)temp[0];
        if (u == NULL) {
          lrwork = 5*MIN(m,n);
        } else {
          /* FIXME: simplify. */
          lrwork = MIN(m,n)*MAX(5*MIN(m,n) + 7, 2*MAX(m,n) + 2*MIN(m,n) + 1);
        }
        offset = ROUND_UP((lrwork + 2*lwork)*sizeof(double), sizeof(INTEGER));
        size = offset + 8*MIN(m,n)*sizeof(INTEGER);
        rwork = ypush_scratch(size, NULL);
        work = INCR_ADDRESS(rwork, lrwork*sizeof(double));
        iwork = INCR_ADDRESS(rwork, offset);
        ZGESDD(jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
               work, &lwork, rwork, iwork, &info);
      }
    }
  } else {
    /* Simple version of the algorithm.*/
    CHARACTER jobu[2], jobvt[2];
    jobu[0] = (u == NULL ? 'N' : (full ? 'A' : 'S'));
    jobu[1] = '\0';
    jobvt[0] = (vt == NULL ? 'N' : (full ? 'A' : 'S'));
    jobvt[1] = '\0';
    if (type == Y_FLOAT) {
      float temp;
      SGESVD(jobu, jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
             &temp, &lwork, &info);
      if (info == 0) {
        lwork = (INTEGER)temp;
        work = ypush_scratch(lwork*sizeof(float), NULL);
        SGESVD(jobu, jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
               work, &lwork, &info);
      }
    } else if (type == Y_DOUBLE) {
      double temp;
      DGESVD(jobu, jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
             &temp, &lwork, &info);
      if (info == 0) {
        lwork = (INTEGER)temp;
        work = ypush_scratch(lwork*sizeof(double), NULL);
        DGESVD(jobu, jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
               work, &lwork, &info);
      }
    } else {
      double temp[2];
      ZGESVD(jobu, jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
             temp, &lwork, NULL, &info);
      if (info == 0) {
        long lrwork = 5*MIN(m,n);
        lwork = (INTEGER)temp[0];
        work = ypush_scratch((2*lwork + lrwork)*sizeof(double), NULL);
        rwork = INCR_ADDRESS(work, 2*lwork*sizeof(double));
        ZGESVD(jobu, jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
               work, &lwork, rwork, &info);
      }
    }
  }

  /* Build result. */
  yarg_drop(1); /* drop workspace */
  if (info == 0) {
    /* Algorithm was successful.  Extract results from the stack. */
    if (vt_ref >= 0L) {
      yput_global(vt_ref, 0);
    }
    if (u_ref >= 0L) {
      yput_global(u_ref, 1);
    }
    if (s_ref >= 0L) {
      yput_global(s_ref, 2);
    }
  } else {
    /* Algorithm failed.  raise an error if called as a subroutine. */
    if (yarg_subroutine()) {
      char *msg;
      if (info > 0) {
        msg = "algorithm did not converge";
      } else {
        sprintf(msgbuf,
                "bug in %s wrapper: %d-th argument had an illegal value",
                (divide_and_conquer ? "xGESDD" : "xGESVD"), -(int)info);
        msg = msgbuf;
      }
      y_error(msg);
    }
  }
  ypush_long(info);
}

typedef void *(push_array_func)(long *dims);

static push_array_func *push_array_table[] = {
  (push_array_func *)ypush_c,
  (push_array_func *)ypush_s,
  (push_array_func *)ypush_i,
  (push_array_func *)ypush_l,
  (push_array_func *)ypush_f,
  (push_array_func *)ypush_d,
  (push_array_func *)ypush_z
};

static size_t elsize_table[] = {
  sizeof(char),
  sizeof(short),
  sizeof(int),
  sizeof(long),
  sizeof(float),
  sizeof(double),
  2*sizeof(double)
};

/* Force a numerical array stored on the stack to be a temporary one. */
static void *force_scratch(int iarg, void *ptr, int type,
                           long ntot, long *dims)
{
  void *tmp;
  if (yarg_scratch(iarg)) {
    return ptr;
  }
  if (type < 0 || type > 6) {
    y_error("bad type");
    return NULL;
  }
  tmp = (*push_array_table[type])(dims);
  memcpy(tmp, ptr, elsize_table[type]*ntot);
  yarg_swap(iarg + 1, 0);
  yarg_drop(1);
  return tmp;
}

static void *push_1d(int type, long n)
{
  long dims[2];
  if (type < 0 || type > 6) {
    y_error("bad type");
    return NULL;
  }
  dims[0] = 1;
  dims[1] = n;
  return (*push_array_table[type])(dims);
}

static void *push_2d(int type, long m, long n)
{
  long dims[3];
  if (type < 0 || type > 6) {
    y_error("bad type");
    return NULL;
  }
  dims[0] = 2;
  dims[1] = m;
  dims[2] = n;
  return (*push_array_table[type])(dims);
}

void Y_lpk_ggsvd(int argc)
{
  array_t a, b;
  long k_ref, l_ref, a_ref, b_ref, alpha_ref, beta_ref;
  long u_ref, v_ref, q_ref, i_ref;
  void *alpha, *beta, *u, *v, *q, *work, *rwork;
  long dims[Y_DIMSIZE], *lwork;
  push_array_t *push_array;
  char jobu[2], jobv[2], jobq[2];
  INTEGER i, m, n, p, k, l, info, ldu, ldv, ldq, *iwork, wn;
  int a_scratch, b_scratch, type, iarg = argc;
  size_t size, rwork_offset, iwork_offset;

  /* Get the arguments. */

  if (--iarg < 0) goto bad_nargs;
  m = ygets_l(iarg);

  if (--iarg < 0) goto bad_nargs;
  n = ygets_l(iarg);

  if (--iarg < 0) goto bad_nargs;
  p = ygets_l(iarg);

  if (--iarg < 0) goto bad_nargs;
  k_ref = yget_ref(iarg);

  if (--iarg < 0) goto bad_nargs;
  l_ref = yget_ref(iarg);

  if (--iarg < 0) goto bad_nargs;
  a_ref = yget_ref(iarg);
  a_scratch = yarg_scratch(iarg);
  get_array(&a, iarg);

  if (--iarg < 0) goto bad_nargs;
  b_ref = yget_ref(iarg);
  b_scratch = yarg_scratch(iarg);
  get_array(&b, iarg);

  if (--iarg < 0) goto bad_nargs;
  alpha_ref = yget_ref(iarg);

  if (--iarg < 0) goto bad_nargs;
  beta_ref = yget_ref(iarg);

#define GET_NEXT_OPTIONAL_REF (--iarg < 0 ? -1L : yget_ref(iarg))
  u_ref = GET_NEXT_OPTIONAL_REF;
  v_ref = GET_NEXT_OPTIONAL_REF;
  q_ref = GET_NEXT_OPTIONAL_REF;
  i_ref = GET_NEXT_OPTIONAL_REF;
#undef GET_NEXT_OPTIONAL_REF

  /* Check arguments. */
  if (m <= 0 || n <= 0 || p <= 0 || m*n != a.ntot || p*n != b.ntot) {
    y_error("bad dimensions M, N and P");
  }
  if (a.type <= Y_FLOAT && b.type <= Y_FLOAT) {
    type = Y_FLOAT;
    push_array = (push_array_t *)ypush_f;
  } else if (a.type == Y_COMPLEX || b.type == Y_COMPLEX) {
    type = Y_COMPLEX;
    push_array = (push_array_t *)ypush_z;
  } else {
    type = Y_DOUBLE;
    push_array = (push_array_t *)ypush_d;
  }

  /* Free memory. */
#define FREE_GLOBAL(ref) if (ref >= -1L) yput_global(ref, 0)
  ypush_nil();
  FREE_GLOBAL(k_ref);
  FREE_GLOBAL(l_ref);
  FREE_GLOBAL(alpha_ref);
  FREE_GLOBAL(beta_ref);
  FREE_GLOBAL(u_ref);
  FREE_GLOBAL(v_ref);
  FREE_GLOBAL(q_ref);
  FREE_GLOBAL(i_ref);
  yarg_drop(1);
#undef FREE_GLOBAL

  /* Convert/copy input arrays. */
  if (a.type != type) {
    coerce_array(&a, type);
  } else if (! a_scratch) {
    copy_array(&a);
  }
  if (b.type != type) {
    coerce_array(&b, type);
  } else if (! b_scratch) {
    copy_array(&b);
  }

  /* Push real arrays ALPHA and BETA. */
  dims[0] = 1;
  dims[1] = n;
  if (type == Y_FLOAT) {
    alpha = ypush_f(dims);
    beta = ypush_f(dims);
  } else {
    alpha = ypush_d(dims);
    beta = ypush_d(dims);
  }

  /* Allocate output arrays. */
  if (u_ref >= 0L) {
    dims[0] = 2;
    dims[1] = m;
    dims[2] = m;
    u = push_array(dims);
    ldu = m;
    jobu[0] = 'U';
  } else {
    u = NULL;
    ldu = 1;
    jobu[0] = 'N';
  }
  jobu[1] = '\0';
  if (v_ref >= 0L) {
    dims[0] = 2;
    dims[1] = p;
    dims[2] = p;
    v = push_array(dims);
    ldv = p;
    jobv[0] = 'V';
  } else {
    v = NULL;
    ldv = 1;
    jobv[0] = 'N';
  }
  jobv[1] = '\0';
  if (q_ref >= 0L) {
    dims[0] = 2;
    dims[1] = n;
    dims[2] = n;
    q = push_array(dims);
    ldq = n;
    jobq[0] = 'Q';
  } else {
    q = NULL;
    ldq = 1;
    jobq[0] = 'N';
  }
  jobq[1] = '\0';

  /* Compute WN = max(3*N,M,P) + N and allocate workspace for WORK and RWORK
     arrays. */
  wn = 3*n;
  if (wn < m) wn = m;
  if (wn < p) wn = p;
  wn += n;
  if (type == Y_COMPLEX) {
    rwork_offset = 2*wn*sizeof(double);
    size = rwork_offset + 2*n*sizeof(double);
  } else {
    rwork_offset = 0;
    if (type == Y_FLOAT) {
      size = wn*sizeof(float);
    } else {
      size = wn*sizeof(double);
    }
  }
  iwork_offset = ROUND_UP(size, sizeof(INTEGER));
  size = iwork_offset + n*sizeof(INTEGER);
  work = ypush_scratch(size, NULL);
  if (type == Y_COMPLEX) {
    rwork = (void *)((char *)work + rwork_offset);
  } else {
    rwork = NULL;
  }
  iwork = (INTEGER *)((char *)work + iwork_offset);

  /* Perform operation. */
  if (type == Y_FLOAT) {
    (void)SGGSVD(jobu, jobv, jobq, &m, &n, &p, &k, &l,
                 a.data, &m, b.data, &p, alpha, beta,
                 u, &ldu, v, &ldv, q, &ldq,
                 work, iwork, &info);
  } else if (type == Y_DOUBLE) {
    (void)DGGSVD(jobu, jobv, jobq, &m, &n, &p, &k, &l,
                 a.data, &m, b.data, &p, alpha, beta,
                 u, &ldu, v, &ldv, q, &ldq,
                 work, iwork, &info);
  } else {
    (void)ZGGSVD(jobu, jobv, jobq, &m, &n, &p, &k, &l,
                 a.data, &m, b.data, &p, alpha, beta,
                 u, &ldu, v, &ldv, q, &ldq,
                 work, rwork, iwork, &info);
  }

  /* Save outputs.  Note that U, V, Q, and I have been pushed on top of the
     stack. */
  if (i_ref >= 0L) {
    dims[0] = 1;
    dims[1] = n;
    lwork = ypush_l(dims);
    for (i = 0; i < n; ++i) {
      lwork[i] = iwork[i];
    }
    yput_global(i_ref, 0);
    yarg_drop(2); /* drop temporary workspace and output LWORK */
  } else {
    yarg_drop(1); /* drop temporary workspace */
  }
#define SAVE_GLOBAL(ref)                        \
  if ((ref) >= 0L) {                            \
    yput_global((ref), 0);                      \
    yarg_drop(1);                               \
  }
  SAVE_GLOBAL(q_ref);
  SAVE_GLOBAL(v_ref);
  SAVE_GLOBAL(u_ref);
  SAVE_GLOBAL(beta_ref);
  SAVE_GLOBAL(alpha_ref);
  SAVE_GLOBAL(b_ref);
  SAVE_GLOBAL(a_ref);
#undef SAVE_GLOBAL
#define SAVE_GLOBAL(ref, val)                   \
  if ((ref) >= 0L) {                            \
    ypush_long(val);                            \
    yput_global((ref), 0);                      \
    yarg_drop(1);                               \
  }
  SAVE_GLOBAL(l_ref, l);
  SAVE_GLOBAL(k_ref, k);
#undef SAVE_GLOBAL

  return;

 bad_nargs:
  y_error("bad number of arguments");
}

/*---------------------------------------------------------------------------*/
/* UTILITIES */

static long get_size_of_linear_system(const long adims[], const long bdims[])
{
  long j, n, s;
  if (adims[0] % 2L != 0L) {
    not_square:
    y_error("left hand side matrix A must be square");
    return -1L;
  }
  s = adims[0]/2L;
  if (bdims[0] < s) {
  incompatible:
    y_error("left hand side matrix A and right hand side vector B have incompatible dimensions");
    return -1L;
  }
  n = 1L;
  for (j = 1; j <= s; ++j) {
    if (adims[j] != bdims[j]) {
      goto incompatible;
    }
    if (adims[j] != adims[s + j]) {
      goto not_square;
    }
    n *= adims[j];
  }
  return n;
}

static long get_linear_system(linear_system_t *sys,
                              int a_iarg, int b_iarg,
                              unsigned int flags)
{
  long a_dims[Y_DIMSIZE], b_dims[Y_DIMSIZE];
  void *a_ptr, *b_ptr;
  long a_ntot, b_ntot, n;
  int type, a_type, b_type, a_temp, b_temp;
  int max_type = ((flags & COMPLEX_SYSTEM) != 0 ? Y_COMPLEX : Y_DOUBLE);

  /* left hand side matrix A */
  a_temp = (yget_ref(a_iarg) == -1L);
  a_ptr = ygeta_any(a_iarg, &a_ntot, a_dims, &a_type);
  if (a_type < Y_CHAR || a_type > max_type) {
    y_error("bad data type for left hand side matrix A");
  }

  /* right hand side vector b */
  b_temp = (yget_ref(b_iarg) == -1L);
  b_ptr = ygeta_any(b_iarg, &b_ntot, b_dims, &b_type);
  if (b_type < Y_CHAR || b_type >  max_type) {
    y_error("bad data type for right hand side vector B");
  }

  /* get the size of the problem */
  n = get_size_of_linear_system(a_dims, b_dims);

  /* convert arrays and prepare result (operation is done in-place) */
  if ((flags & REAL_OR_COMPLEX_SYSTEM) == COMPLEX_SYSTEM) {
    type = Y_COMPLEX;
  } else {
    type = MAX(a_type, b_type);
    if (type < Y_FLOAT && type >= Y_CHAR) {
      type = Y_FLOAT;
    }
  }
  if (a_type != type) {
    a_ptr = ygeta_coerce(a_iarg, a_ptr, a_ntot, a_dims, a_type, type);
    a_temp = 1;
  } else if (! a_temp) {
    if (type == Y_FLOAT) {
      a_ptr = push_copy_f(a_ptr, a_dims);
    } else if (type == Y_DOUBLE) {
      a_ptr = push_copy_d(a_ptr, a_dims);
    } else /* Y_COMPLEX */ {
      a_ptr = push_copy_z(a_ptr, a_dims);
    }
    yarg_swap(0, a_iarg + 1);
    yarg_drop(1);
  }
  if (b_type != type) {
    b_ptr = ygeta_coerce(b_iarg, b_ptr, b_ntot, b_dims, b_type, type);
    b_temp = 1;
  } else if (! b_temp) {
    if (type == Y_FLOAT) {
      b_ptr = push_copy_f(b_ptr, b_dims);
    } else if (type == Y_DOUBLE) {
      b_ptr = push_copy_d(b_ptr, b_dims);
    } else /* Y_COMPLEX */ {
      b_ptr = push_copy_z(b_ptr, b_dims);
    }
    yarg_swap(0, b_iarg + 1);
    yarg_drop(1);
  }

  sys->a = a_ptr;
  sys->b = b_ptr;
  sys->type = type;
  sys->lda = sys->ldb = sys->n = n;
  sys->nrhs = b_ntot/n;
  return n;
}

static float get_f(int iarg)
{
  int type;
  long dims[Y_DIMSIZE];
  void *p = ygeta_any(iarg, NULL, dims, &type);
  if (dims[0] == 0) {
    switch (type) {
#define CASE(id,ctype) case id: return (float)*(ctype *)p
      CASE(Y_CHAR,   char);
      CASE(Y_SHORT,  short);
      CASE(Y_INT,    int);
      CASE(Y_LONG,   long);
      CASE(Y_FLOAT,  float);
      CASE(Y_DOUBLE, double);
#undef CASE
    }
  }
  y_error("expecting a non-complex numerical scalar");
  return 0.0f;
}

static double get_d(int iarg)
{
  int type;
  long dims[Y_DIMSIZE];
  void *p = ygeta_any(iarg, NULL, dims, &type);
  if (dims[0] == 0) {
    switch (type) {
#define CASE(id,ctype) case id: return (double)*(ctype *)p
      CASE(Y_CHAR,   char);
      CASE(Y_SHORT,  short);
      CASE(Y_INT,    int);
      CASE(Y_LONG,   long);
      CASE(Y_FLOAT,  float);
      CASE(Y_DOUBLE, double);
#undef CASE
    }
  }
  y_error("expecting a non-complex numerical scalar");
  return 0.0;
}

static void get_z(int iarg, double z[2])
{
  int type;
  long dims[Y_DIMSIZE];
  void *p = ygeta_any(iarg, NULL, dims, &type);
  if (dims[0] == 0) {
    switch (type) {
#define CASE(id,ctype) case id: z[0] = *(ctype *)p; z[1] = 0.0; return
      CASE(Y_CHAR,   char);
      CASE(Y_SHORT,  short);
      CASE(Y_INT,    int);
      CASE(Y_LONG,   long);
      CASE(Y_FLOAT,  float);
      CASE(Y_DOUBLE, double);
#undef CASE
    case Y_COMPLEX:
      z[0] = ((double *)p)[0];
      z[1] = ((double *)p)[1];
      return;
    }
  }
  y_error("expecting a complex scalar");
}

#ifdef USE_CBLAS

static CBLAS_DIAG get_cblas_diag(int iarg)
{
  switch (ygets_i(iarg)) {
  case LPK_NON_UNIT: return CblasNonUnit;
  case LPK_UNIT:     return CblasUnit;
  }
  y_error("invalid value for DIAG");
  return CblasNonUnit; /* to avoid compiler warnings */
}

static CBLAS_ORDER get_cblas_order(int iarg)
{
  switch (ygets_i(iarg)) {
  case LPK_ROW_MAJOR: return CblasRowMajor;
  case LPK_COL_MAJOR: return CblasColMajor;
  }
  y_error("invalid value for ORDER");
  return CblasColMajor; /* to avoid compiler warnings */
}

static CBLAS_SIDE get_cblas_side(int iarg)
{
  switch (ygets_i(iarg)) {
  case LPK_LEFT:  return CblasLeft;
  case LPK_RIGHT: return CblasRight;
  }
  y_error("invalid value for SIDE");
  return CblasLeft; /* to avoid compiler warnings */
}

static CBLAS_TRANSPOSE get_cblas_trans(int iarg)
{
  switch (ygets_i(iarg)) {
  case LPK_NO_TRANS:   return CblasNoTrans;
  case LPK_TRANS:      return CblasTrans;
  case LPK_CONJ_TRANS: return CblasConjTrans;
  }
  y_error("invalid value for TRANS");
  return CblasNoTrans; /* to avoid compiler warnings */
}

static CBLAS_UPLO get_cblas_uplo(int iarg)
{
  switch (ygets_i(iarg)) {
  case LPK_UPPER: return CblasUpper;
  case LPK_LOWER: return CblasLower;
  }
  y_error("invalid value for UPLO");
  return CblasUpper; /* to avoid compiler warnings */
}

#endif /* USE_CBLAS */

static void get_lapack_trans(int iarg, CHARACTER trans[2])
{
  switch (ygets_i(iarg)) {
  case LPK_NO_TRANS:   trans[0] = 'N'; break;
  case LPK_TRANS:      trans[0] = 'T'; break;
  case LPK_CONJ_TRANS: trans[0] = 'C'; break;
  default: y_error("invalid value for TRANS");
  }
  trans[1] = '\0';
}

static void get_lapack_uplo(int iarg, CHARACTER uplo[2])
{
  switch (ygets_i(iarg)) {
  case LPK_UPPER: uplo[0] = 'U'; break;
  case LPK_LOWER: uplo[0] = 'L'; break;
  default: y_error("invalid value for UPLO");
  }
  uplo[1] = '\0';
}

static void push_string(const char *value)
{
  ypush_q((long *)NULL)[0] = (value ? p_strcpy((char *)value) : NULL);
}

static float *push_copy_f(const float *src, const long dims[])
{
  float *dst = ypush_f((long *)dims);
  long i, ntot = 1;
  for (i = dims[0]; i > 0; --i) {
    ntot *= dims[i];
  }
  return memcpy(dst, src, ntot*sizeof(float));
}

static double *push_copy_d(const double *src, const long dims[])
{
  double *dst = ypush_d((long *)dims);
  long i, ntot = 1;
  for (i = dims[0]; i > 0; --i) {
    ntot *= dims[i];
  }
  return memcpy(dst, src, ntot*sizeof(double));
}

static double *push_copy_z(const double *src, const long dims[])
{
  double *dst = ypush_z((long *)dims);
  long i, ntot = 1;
  for (i = dims[0]; i > 0; --i) {
    ntot *= dims[i];
  }
  return memcpy(dst, src, 2*ntot*sizeof(double));
}

static void copy_array(array_t *arr)
{
  long elem_size;
  void *new_data;
  switch (arr->type) {
  case Y_CHAR:
    new_data = ypush_c(arr->dims);
    elem_size = sizeof(char);
    break;
  case Y_SHORT:
    new_data = ypush_s(arr->dims);
    elem_size = sizeof(short);
    break;
  case Y_INT:
    new_data = ypush_i(arr->dims);
    elem_size = sizeof(int);
    break;
  case Y_LONG:
    new_data = ypush_l(arr->dims);
    elem_size = sizeof(long);
    break;
  case Y_FLOAT:
    new_data = ypush_f(arr->dims);
    elem_size = sizeof(float);
    break;
  case Y_DOUBLE:
    new_data = ypush_d(arr->dims);
    elem_size = sizeof(double);
    break;
  case Y_COMPLEX:
    new_data = ypush_d(arr->dims);
    elem_size = 2*sizeof(double);
    break;
  default:
    return;
  }
  arr->data = memcpy(new_data, arr->data, arr->ntot*elem_size);
  yarg_swap(arr->iarg + 1, 0);
  yarg_drop(1);
}

static void get_array(array_t *arr, int iarg)
{
  arr->data = ygeta_any(iarg, &arr->ntot, arr->dims, &arr->type);
  if (arr->type > Y_COMPLEX) {
    y_error("expecting numerical array");
    return;
  }
  arr->iarg = iarg;
}

static void get_vector(vector_t *vec, int iarg, int rng)
{
  vec->arr.data = ygeta_any(iarg, &vec->arr.ntot, vec->arr.dims,
                            &vec->arr.type);
  if (vec->arr.type > Y_COMPLEX) {
    y_error("expecting numerical array");
    return;
  }
  vec->arr.iarg = iarg;
  if (rng < 0) {
    vec->offset = 0;
    vec->size = vec->arr.ntot;
    vec->step = 1;
    vec->rank = vec->arr.dims[0];
  } else {
    long range[3];
    long start, stop, step, size, length = vec->arr.ntot;
    int flags = yget_range(rng, range);
    if (flags == 0 || (flags & (1|Y_MIN_DFLT|Y_MAX_DFLT)) != flags) {
      y_error("expecting a range in the form START:STOP:STEP");
      return;
    }
    step = range[2];
#if 0 /* STEP = 0 is forbidden in Yorick */
    if (step == 0) {
      y_error("0 step in index range (start:stop:step)");
      return;
    }
#endif
    if ((flags & Y_MIN_DFLT) == 0) {
      start = range[0];
      if (start <= 0) start += length; /* apply Yorick's rule */
      if (start < 1 || start > length) {
        y_error("out of range start index");
        return;
      }
    } else {
      start = (step < 0 ? length : 1);
    }
    if ((flags & Y_MAX_DFLT) == 0) {
      stop = range[1];
      if (stop <= 0) stop += length; /* apply Yorick's rule */
      if (stop < 1 || stop > length) {
        y_error("out of range stop index");
        return;
      }
    } else {
      stop = (step < 0 ? 1 : length);
    }
    if (step < 0) {
      if (start < stop) {
        goto wrongSign;
      }
      size = (start - stop)/(-step) + 1L;
      vec->offset = start - 1L + (size - 1L)*step;
    } else {
      if (start > stop) {
        goto wrongSign;
      }
      size =  (stop - start)/step + 1L;
      vec->offset = start - 1L;
    }
    vec->step = step;
    vec->size = size;
    vec->rank = 1;
  }
  return;
 wrongSign:
  y_error("array index range step has wrong sign");
}

static void get_matrix(matrix_t *mat, int iarg)
{
  long rank;
  mat->arr.data = ygeta_any(iarg, &mat->arr.ntot, mat->arr.dims,
                            &mat->arr.type);
  if (mat->arr.type > Y_COMPLEX) {
    y_error("expecting numerical array");
    return;
  }
  mat->arr.iarg = iarg;
  rank = mat->arr.dims[0];
  mat->nrows = (rank >= 1L ? mat->arr.dims[1] : 1);
  mat->ncols = (rank >= 2L ? mat->arr.ntot/mat->arr.dims[1] : 1);
}

static BLAS_INTEGER check_vector_sizes(const vector_t *x, const vector_t *y)
{
  int i, rank;
  if ((rank = x->rank) != y->rank || x->size != y->size) {
  badSize:
    y_error("incompatible vector dimensions or ranges");
    return -1;
  }
  if (rank > 1) { /* no needs to check dimension lists for flat vectors */
    for (i = 1; i <= rank; ++i) {
      if (x->arr.dims[i] != y->arr.dims[i]) {
        goto badSize;
      }
    }
  }
  return x->size;
}

static void coerce_array(array_t *arr, int type)
{
  arr->data = ygeta_coerce(arr->iarg, arr->data, arr->ntot,
                           arr->dims, arr->type, type);
  arr->type = type;
}

static void coerce_vector(vector_t *vec, int type)
{
  vec->arr.data = ygeta_coerce(vec->arr.iarg, vec->arr.data, vec->arr.ntot,
                               vec->arr.dims, vec->arr.type, type);
  vec->arr.type = type;
}

static void coerce_matrix(matrix_t *mat, int type)
{
  mat->arr.data = ygeta_coerce(mat->arr.iarg, mat->arr.data, mat->arr.ntot,
                               mat->arr.dims, mat->arr.type, type);
  mat->arr.type = type;
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * c-basic-offset: 2
 * fill-column: 78
 * coding: utf-8
 * End:
 */
