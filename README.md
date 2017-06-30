               __   ___                            _
               \ \ / / |     __ _ _ __   __ _  ___| | __
                \ V /| |    / _` | '_ \ / _` |/ __| |/ /
                 | | | |___| (_| | |_) | (_| | (__|   <
                 |_| |_____|\__,_| .__/ \__,_|\___|_|\_\
                                 |_|

YLapack is a semi-low level Yorick interface to BLAS and LAPACK libraries.
The reasons for this interface are multiples:

1. Yorick only provides general matrix inversion and resolution of linear
   systems of equations (LUsolve), resolution of tridiagonal systems of
   equations (TDsolve), solving linear least squares (QRsolve or SVsolve) and
   singular value deconposition (SVdec).  Lapack gives many more: eigenvalues
   and eigenvectors, symmetric matrices, etc.

2. There are a number of faster implementations of BLAS and LAPACK than the one
   built into Yorick (an automatic Fortran to C conversion by a smart
   Emacs-lisp script written by Dave Munro plus hand editing).  You may also
   want to benefit from the speed-up of basic linear algebra operations such as
   applying a matrix-vector or matrix-matrix multiplication, computing the
   scalar product of two vectors (sum(x*y) in Yorick), computing the Euclidean
   norm of a vector, etc.

3. Although it should be possible to link Yorick with these fast libraries to
   benefit from the speed in LUsolve, QRsolve, TDsolve, SVsolve and SVdec, you
   may also want to have a fancier interface to these functions which
   generalizes the rules for the dot product so that it can be applied on
   several consecutive dimensions at the same time.

The current version of YLapack is merely a proof of concept to check whether
Yorick can call multi-threaded functions and benefit from the speed up of one
of the LAPACK implementation for your machine: GotoBlas, MKL, Atlas, etc.  The
main limitation is that only a subset of BLAS and LAPACK functions have been
interfaced.


## PORTABILITY

There are a lot of macro definitions to make this software as portable as
possible, to link with different implementations of LAPACK and even compile
the plugin with a different compiler than the one used for Yorick itself.  For
instance, I was able to run the tests with a plugin compiled with ICC (the
Intel C compiler) and linked with the MKL inside a Yorick compiled with GCC
(the GNU C compiler).


## DESIGN

LAPACK functions which return a status (generally in an integer variable
called INFO) are mapped to a Yorick function which returns the status, it is
the caller's responsibility to check the returned value (non-zero means error)
however when called as a subroutine, as there are no means for the caller to
realize that an error occured, the Yorick wrapper will raise an error.

Dimensions are automatically guessed from the arguments.


## LINKING WITH GOTOBLAS2

[GotoBLAS2](http://www.tacc.utexas.edu/tacc-projects/gotoblas2) is a very fast
BLAS library (plus some LAPACK functions) developped by Kazushige Goto for a
number of processors.  You can download the last release made by Kazushige
Goto (BSD license) at:

    http://cms.tacc.utexas.edu/fileadmin/images/GotoBLAS2-1.13_bsd.tar.gz

Unfortunately, GotoBLAS2 is no longer maintained by its author.  To my
knowledge, there are two open-source projects which aim at maintaining
GotoBLAS2:

* OpenBLAS at http://xianyi.github.com/OpenBLAS/

* SurviveGotoBLAS2 at http://prs.ism.ac.jp/~nakama/SurviveGotoBLAS2/

They are worth having a look, especially if you have some recent processor not
supported by GotoBLAS2-1.13.  They also fix a number of bugs in GotoBLAS2.  We
currently use the OpenBLAS (0.2.5) variant with no problem.  See notes below if
you want to use the SurviveGotoBLAS2 variant.

After unpacking the archive, enter the source directory of GotoBLAS2 and just
type:

    shell> make

to build the library for your processor and compiler; if you want to support
multiple architecture, build the library with:

    shell> make DYNAMIC_ARCH=1

If you want the dynamic library:

    shell> make [OPTIONS] shared

where `OPTIONS` are the same (e.g., `DYNAMIC_ARCH=1`) as used to build the
static library.  Note that the static library is build as position independent
code (PIC), so it is perfectly usable for making a Yorick plugin.

If you plan to use multi-threaded FFTW, add the options:

    USE_THREAD=1 USE_OPENMP=1

when building GotoBLAS2.

Note that to successfully compile for multiple architectures the
GotoBLAS2-1.13_bsd version, I had to patch the file `driver/others/dynamic.c`,
the differences are:

```
71c71
< static int get_vendor(void){
---
> static int get_vendor(void) {
77,79c77,79
<   *(int *)(&vendor[0]) = ebx;
<   *(int *)(&vendor[4]) = edx;
<   *(int *)(&vendor[8]) = ecx;
---
>   memcpy(&vendor[0], &ebx, 4);
>   memcpy(&vendor[4], &edx, 4);
>   memcpy(&vendor[8], &ecx, 4);
201c201
<   if (gotoblas == NULL) gotoblas = gotoblas_KATMAI;
---
>   if (gotoblas == NULL) gotoblas = &gotoblas_KATMAI;
203c203
<   if (gotoblas == NULL) gotoblas = gotoblas_PRESCOTT;
---
>   if (gotoblas == NULL) gotoblas = &gotoblas_PRESCOTT;
```

These have been fixed in SurviveGotoBLAS2 which has the advantage of taking
into account newest LAPACK version (3.3.1 as of writing) while
GotoBLAS2-1.13_bsd is stuck at LAPACK version 3.1.1.

To build SurviveGotoBLAS2, download the latest version, unpack it, and:

    shell> make DYNAMIC_ARCH=1 LAPACK_VERSION=3.3.1

I do not use options `NO_WARMUP=1 NO_AFFINITY=1 NUM_THREADS=48` since I got a
segmentation fault in the plugin (in lpk_gesv) when using the library compiled
with these options (though the tests were successfully passed and I do not
know which of this option is responsible of the problem).  Avoid
DYNAMIC_ARCH=1 if you just want a version for your machine.

Once you have built the GotoBLAS library, copy it (or make a symbolic link or
change the rules in `rules/Make.gotoblas2`) into the directory of YLapack
source as `libgoto2.a` and run:

    shell> yorick -batch make.i
    shell> make clean
    shell> make MODEL=gotoblas2

then you can test the plugin (if you have used MKL or shared libraries, you
may have to set `LD_LIBARY_PATH` accordingly):

    shell> yorick

    yorick> include, "lapack-test.i";
    yorick> lpk_test_gesv, 3000, nloops=20;

and enjoy the speedup ;-)


## PERFORMANCES

Execution times on my laptop -- Intel Core i7 (Q820) at 1.73GHz -- (the
percentage is the rate of CPU occupation, more than 100% means
multi-threading; the speed-up w.r.t. Yorick is given between the square
brackets).

### Scalar product
```
-----------------------------------------------------------------------
             Yorick            GotoBlas2              MKL
   size     sum(x*y)           lpk_dot              lpk_dot
                          (µs = microseconds)
-----------------------------------------------------------------------
   10,000    21 µs (100%)   5.4 µs (100%) [4.0]     6.7 µs (396%) [3.1]
  100,000   220 µs (100%)    55 µs (100%) [4.1]      46 µs (393%) [4.8]
1,000,000  3700 µs (100%)  1200 µs (100%) [3.1]     948 µs (398%) [3.9]
  1024^2   4000 µs (100%)  1300 µs (100%) [3.1]    1000 µs (399%) [4.0]
-----------------------------------------------------------------------
```
For this BLAS level 1 operation, GotoBlas2 is not multi-threaded.
MKL is the fastest (for vectors of size > 100,000).  MKL and GotoBlas2
provide speed-up between 3 and 5.


### Resolution of a linear system of equations
```
----------------------------------------------------------------------
            Yorick            GotoBlas2              MKL
 size      LUsolve            lpk_gesv             lpk_gesv
----------------------------------------------------------------------
   100  0.39 ms  (99%)   0.25 ms (189%)  [1.6]   0.22 ms (392%)  [1.7]
   500    39 ms  (99%)     10 ms (212%)  [3.9]    6.1 ms (397%)  [6.3]
 1,000   0.30 s (100%)   0.066 s (338%)  [4.5]   0.038 s (397%)  [8.1]
 2,000    2.3 s (100%)    0.27 s (355%)  [8.6]    0.27 s (398%)  [8.6]
 3,000    8.0 s (100%)    0.86 s (355%)  [9.4]    0.91 s (387%)  [8.9]
 5,000     38 s (100%)     3.8 s (367%) [10.1]     4.0 s (392%)  [9.6]
10,000    303 s (100%)      30 s (379%) [10.2]      28 s (392%) [11.1]
----------------------------------------------------------------------
```
Hence for moderate size matrix MKL is the fastest but for very large
matrices GotoBlas2 and MKL have similar speed-up of ~ 10.

### Tests with GotoBLAS2 on Linux with CPU Intel Core i7-2600 @ 3.40GHz
```
---------------------------------------------------------------------
                  Yorick                  Lapack
 size             LUsolve                lpk_gesv
---------------------------------------------------------------------
     50  3.53E-05 +/- 6E-06 ( 96%)   2.99E-05 +/- 4E-05 (139%)  [1.2]
    100  2.41E-04 +/- 1E-05 (100%)   9.75E-05 +/- 5E-06 (156%)  [2.5]
    200  1.76E-03 +/- 3E-05 ( 99%)   4.37E-04 +/- 2E-05 (178%)  [4.0]
    300  5.51E-03 +/- 8E-05 ( 99%)   1.12E-03 +/- 2E-04 (204%)  [4.9]
    500  2.41E-02 +/- 4E-04 (100%)   4.30E-03 +/- 1E-04 (216%)  [5.6]
  1,000  1.90E-01 +/- 9E-04 (100%)   2.12E-02 +/- 4E-04 (305%)  [9.0]
  2,000  1.47E+00 +/- 1E-03 (100%)   1.33E-01 +/- 2E-03 (358%) [11.0]
  3,000  5.01E+00 +/- 4E-03 ( 99%)   4.34E-01 +/- 5E-03 (353%) [11.5]
  5,000  2.32E+01 +/- 1E-02 ( 99%)   1.83E+00 +/- 3E-03 (369%) [12.7]
 10,000  1.84E+02 +/- 6E-02 (100%)   1.36E+01 +/- 7E-02 (382%) [13.5]
 20,000  1.76E+03 +/- 1E+01 ( 99%)   1.04E+02 +/- 8E-02 (390%) [17.0]
---------------------------------------------------------------------
```
Note: I have tested Yorick built-in LUsolve (compiled with GCC 4.5.2)
versus DGESV in Lapack 3.3.1 (compiled with Intel ifort XE 2011 SP1.6.233)
and noticed a speed-up of ~ 1.5 for DGESV.


## Cholesky factorization

Perform Cholesky factorization (DPOTRF) and solve a linear system with a
positive definite symmetric left hand side matrix with GotoBLAS2.
```
----------------------------------------------------------------
         Size                    laptop           server
----------------------------------------------------------------
DPOTRF  5,000x5,000      2.2 sec (320%)     1.3 sec (330%)
DPOTRF 10,000x10,000      (not done)       10.5 sec (320%)
DPOTRF 12,000x12,000     23 sec (387%)     15.6 sec (368%)
----------------------------------------------------------------
DPOTRS 5,000x5,000      0.2 sec
----------------------------------------------------------------
```
* laptop = Linux laptop with Intel Core i7 Q820 at 1.73GHz
* server = Linux server with Intel Xeon X3450 at 2.67Ghz
* dpotrf = Cholesky decomposition (done once for a given C)
* dpotrs = solve the system given the Cholesky decomposition

These times can be compared to LUsolve.


# WORK IN PROGRESS

## LINKING WITH MKL

When linking with the Math Kernel Library (MKL), you need to use 3 or 4
libraries:

1. Interface layer.

   For the IA-32 and Intel(R) MIC targets, there are only one MKL interface
   with 32-bits integers.  For the Intel(R) 64 target, there are two possible
   MKL interfaces: lp64 with 32-bits integers and ilp64 with 64-bits integers.

   - `libmkl_gf_ilp64.a` for GNU Fortran compiler with 64-bit integers on
     64-bit processor

   - `libmkl_intel_ilp64.a` for Intel compiler

2. Threading layer: choose a multi-thread library.

3. Computational layer: a multi-thread library.

4. Run-time libraries (only with MPI?).

To figure out which libraries to use with MKL, you can have a look at "Intel(c)
Math Kernel Library Link Line Advisor":

http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/

To start Yorick with a given `LD_LIBRARY_PATH`:

    LD_LIBRARY_PATH=/opt/intel/lib/intel64:/usr/local/lib rlwrap yorick

