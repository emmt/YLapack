/*
 * ylapack.h --
 *
 * Definitions for BLAS, CBLAS, GOTOBLAS, OPENBLAS and LAPACK.
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

#ifndef _YLAPACK_H
#define _YLAPACK_H 1

/* The following macros should be set to match the data types for the Fortran
   compiler used to build BLAS and/or LAPACK libraries. */
#ifndef CONST
# define CONST     const
#endif
#ifndef CHARACTER
# define CHARACTER char
#endif
#ifndef LOGICAL
# define LOGICAL   int
#endif
#ifndef INTEGER
# if defined(USE_MKL) && defined(MKL_ILP64)
#  define INTEGER   long
# else
#  define INTEGER   int
# endif
#endif
#ifndef VOID
# define VOID      int /* type returned by FORTRAN subroutines */
#endif
#ifndef COMPLEX8
struct _ylpk_complex8 { float z[2]; };
typedef struct _ylpk_complex8 ylpk_complex8;
# define COMPLEX8   ylpk_complex8
#endif
#ifndef COMPLEX16
struct _ylpk_complex16 { double z[2]; };
typedef struct _ylpk_complex16 ylpk_complex16;
# define COMPLEX16  ylpk_complex16
#endif

#ifndef COMPLEX8_PTR
# define COMPLEX8_PTR void *
#endif

#ifndef COMPLEX16_PTR
# define COMPLEX16_PTR void *
#endif

#if defined(USE_GOTOBLAS) || defined(USE_OPENBLAS) || defined(USE_MKL)
# define USE_CBLAS /* GotoBlas and MKL includes a version of CBLAS */
#endif

#if defined(USE_CBLAS) && ! defined(ORDER)
# define ORDER     CblasColMajor /* Yorick is column-major order */
#endif

/* FIXME: BLASINT must be int (32-bit version of the library) or long (64-bit
          version of the library). */
#ifdef USE_64BIT_INT
# define blasint long
#else
# define blasint int
#endif


/*
 * Fortran compilers mangle names differently the name of a FORTRAN symbol may
 * be in upper or lower case letters and may be prefixed or post fixed by an
 * underscore letter.  To cope with this variety, we define the macro
 * FORTRAN_NAME to build the name of a FORTRAN symbol for the C compiler.
 * This macros take two arguments: the FORTRAN symbol name in lowercase and in
 * uppercase letters.  The value of the macro FORTRAN_STYLE may be set to
 * select the appropriate macro definition.  The value of FORTRAN_STYLE is an
 * integer value in the range 0:7 as follows:
 *
 *   - 1st bit set to use uppercase letters, otherwise lowercase;
 *   - 2nd bit set to use an underscore prefix;
 *   - 3rd bit set to use an underscore suffix.
 */
#ifndef FORTRAN_STYLE
# ifdef _WIN32
#  define FORTRAN_STYLE 1
# else
#  define FORTRAN_STYLE 2
# endif
#endif
#if (FORTRAN_STYLE == 0)
# define FORTRAN_NAME(name,NAME) name
#elif (FORTRAN_STYLE == 1)
# define FORTRAN_NAME(name,NAME) NAME
#elif (FORTRAN_STYLE == 2)
# define FORTRAN_NAME(name,NAME) name##_
#elif (FORTRAN_STYLE == 3)
# define FORTRAN_NAME(name,NAME) NAME##_
#elif (FORTRAN_STYLE == 4)
# define FORTRAN_NAME(name,NAME) _##name
#elif (FORTRAN_STYLE == 5)
# define FORTRAN_NAME(name,NAME) _##NAME
#elif (FORTRAN_STYLE == 6)
# define FORTRAN_NAME(name,NAME) _##name##_
#elif (FORTRAN_STYLE == 7)
# define FORTRAN_NAME(name,NAME) _##NAME##_
#else
# error illegal FORTRAN_STYLE
#endif

#undef _YLPK_BEGIN_DECLS
#undef _YLPK_END_DECLS
#ifdef __cplusplus
# define _YLPK_BEGIN_DECLS extern "C" {
# define _YLPK_END_DECLS }
#else
# define _YLPK_BEGIN_DECLS           /* empty */
# define _YLPK_END_DECLS             /* empty */
#endif

#define _YLPK_JOIN2(a1,a2)       a1##a2
#define _YLPK_JOIN3(a1,a2,a3)    a1##a2##a3
#define YLPK_JOIN2(a1,a2)       _YLPK_JOIN2(a1,a2)
#define YLPK_JOIN3(a1,a2,a3)    _YLPK_JOIN3(a1,a2,a3)

/*
 * p - type prefix (s for float, d for double, c for single precision complex,
 *     z for double precision complex, ...)
 * f - function name
 */
#define LAPACK_NAME(name,NAME)  FORTRAN_NAME(name,NAME)
#define CBLAS_NAME(name)        YLPK_JOIN2(cblas_,name)
#define BLAS_NAME(name,NAME)    LAPACK_NAME(name,NAME)

#if defined(USE_MKL)
# include "mkl.h"
#elif defined(USE_CBLAS)
# include "ylapack_cblas.h"
#else
# include "ylapack_blas.h"
#endif
#ifndef USE_MKL
# include "ylapack_lapack.h"
#endif


#endif /* _YLAPACK_H */
