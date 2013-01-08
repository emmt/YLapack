# these values filled in by    yorick -batch make.i
Y_MAKEDIR=/home/apps/libexec/yorick/current
Y_EXE=/home/apps/libexec/yorick/current/bin/yorick
Y_EXE_PKGS=
Y_EXE_HOME=/home/apps/libexec/yorick/current
Y_EXE_SITE=/home/apps/libexec/yorick/current
Y_HOME_PKG=

# ----------------------------------------------------- optimization flags

# options for make command line, e.g.-   make COPT=-g TGT=exe
COPT=$(COPT_DEFAULT)
TGT=$(DEFAULT_TGT)

# ------------------------------------------------ macros for this package

PKG_NAME=ylapack
PKG_I=lapack.i

OBJS = ylapack.o

# Makefile variables LAPACK_DEFS and LAPACK_LIBS must be defined to compile
# and link with your preferred version of LAPACK.
#
#   LAPACK_DEFS is for pre-processor directives
#
#   LAPACK_LIBS is for libraries
#
# in LAPACK_DEFS, you may define some pre-processor macros:
#
#   -DUSE_GOTOBLAS  - to use GotoBlas;
#   -DUSE_MKL       - to use the "Math Kernel Library";
#   -DUSE_CBLAS     - to use C-BLAS interface GotoBlas (set by default when
#                     USE_GOTOBLAS or USE_MKL are defined);
#   -DINTEGER=long  - if a Fortran INTEGER is a "long int" rather than a
#                     just an "int" as assumed by default;
#
MODEL = default

# To use MKL for ia32 processor:
# MKL_DIR = /opt/intel/mkl
# MKL_ARCH = ia32
# MKL_INTERFACE = mkl_gf
# MKL_THREAD = mkl_gnu_thread
# MKL_CORE = mkl_core
# MKL_DEFS = -DUSE_MKL -DINTEGER=int -I$(MKL_DIR)/include
# LAPACK_DEFS = $(MKL_DEFS)
# LAPACK_LIBS = -L$(MKL_DIR)/lib/$(MKL_ARCH) -l$(MKL_INTERFACE) \
#    -l$(MKL_THREAD) -l$(MKL_CORE)

# To use MKL for intel64 processor with 32-bit integers:
# MKL_DIR = /opt/intel/mkl
# MKL_ARCH = intel64
# MKL_INTERFACE = mkl_gf_lp64
# MKL_THREAD = mkl_gnu_thread
# MKL_CORE = mkl_core
# MKL_DEFS = -DUSE_MKL -DINTEGER=int -I$(MKL_DIR)/include
# LAPACK_DEFS = $(MKL_DEFS)
# LAPACK_LIBS = -L$(MKL_DIR)/lib/$(MKL_ARCH) -l$(MKL_INTERFACE) \
#    -l$(MKL_THREAD) -l$(MKL_CORE)

# To use MKL for intel64 processor with 64-bit integers:
# MKL_DIR = /opt/intel/mkl
# MKL_ARCH = intel64
# MKL_INTERFACE = mkl_gf_ilp64
# MKL_THREAD = mkl_gnu_thread
# MKL_CORE = mkl_core
# MKL_DEFS = -DUSE_MKL -DMKL_ILP64 -DINTEGER=long -I$(MKL_DIR)/include
# LAPACK_DEFS = $(MKL_DEFS)
# LAPACK_LIBS = -L$(MKL_DIR)/lib/$(MKL_ARCH) -l$(MKL_INTERFACE) \
#    -l$(MKL_THREAD) -l$(MKL_CORE)


YLAPACK_HEADERS = ylapack.h ylapack_lapack.h ylapack_blas.h ylapack_cblas.h
YLAPACK_RULES = rules/Make.atlas \
                rules/Make.default \
                rules/Make.gotoblas2 \
                rules/Make.icc \
                rules/Make.mkl \
                rules/Make.mkl_dynamic \
                rules/Make.mkl_gnu_ia32 \
                rules/Make.mkl_gnu_ilp64 \
                rules/Make.mkl_gnu_lp64 \
                rules/Make.mkl_icc_ia32 \
                rules/Make.mkl_icc_ilp64 \
                rules/Make.mkl_icc_lp64 \
                rules/Make.mkl_static

# change to give the executable a name other than yorick
PKG_EXENAME=yorick

# PKG_DEPLIBS=-Lsomedir -lsomelib   for dependencies of this package
PKG_DEPLIBS=$(LAPACK_LIBS)
# set compiler (or rarely loader) flags specific to this package
PKG_CFLAGS=$(LAPACK_DEFS)
PKG_LDFLAGS=

# list of additional package names you want in PKG_EXENAME
# (typically Y_EXE_PKGS should be first here)
EXTRA_PKGS=$(Y_EXE_PKGS)

# list of additional files for clean
PKG_CLEAN=

# autoload file for this package, if any
PKG_I_START=
# non-pkg.i include files for this package, if any
PKG_I_EXTRA=

RELEASE_FILES = AUTHORS LICENSE Makefile NEWS README TODO \
	$(PKG_I) $(PKG_I_EXTRA) $(YLAPACK_HEADERS) ylapack.c \
	lapack-test.i rules
RELEASE_NAME = $(PKG_NAME)-$(RELEASE_VERSION).tar.bz2

# -------------------------------- standard targets and rules (in Makepkg)

# set macros Makepkg uses in target and dependency names
# DLL_TARGETS, LIB_TARGETS, EXE_TARGETS
# are any additional targets (defined below) prerequisite to
# the plugin library, archive library, and executable, respectively
PKG_I_DEPS=$(PKG_I)
Y_DISTMAKE=distmake

include $(Y_MAKEDIR)/Make.cfg
include $(Y_MAKEDIR)/Makepkg
include $(Y_MAKEDIR)/Make$(TGT)

include ./rules/Make.$(MODEL)

# override macros Makepkg sets for rules and other macros
# Y_HOME and Y_SITE in Make.cfg may not be correct (e.g.- relocatable)
Y_HOME=$(Y_EXE_HOME)
Y_SITE=$(Y_EXE_SITE)

# reduce chance of yorick-1.5 corrupting this Makefile
MAKE_TEMPLATE = protect-against-1.5

# ------------------------------------- targets and rules for this package

# simple example:
#myfunc.o: myapi.h
# more complex example (also consider using PKG_CFLAGS above):
#myfunc.o: myapi.h myfunc.c
#	$(CC) $(CPPFLAGS) $(CFLAGS) -DMY_SWITCH -o $@ -c myfunc.c

FETCH_DEFS=tclsh ./fetch.tcl
LAPACK_LIBRARY=/usr/lib/atlas-base/atlas/liblapack.a
BLAS_LIBRARY=/usr/lib/atlas-base/atlas/libblas.a
ylapack_lapack.def: fetch.tcl
	$(FETCH_DEFS) --lapack "$(LAPACK_LIBRARY)" > "$@"
ylapack_blas.def: fetch.tcl
	$(FETCH_DEFS) --lapack "$(BLAS_LIBRARY)" > "$@"
ylapack.o: ylapack.c  $(YLAPACK_HEADERS)

test_version.o: test_version.c $(YLAPACK_HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o "$@" -c "$<"
test_version$(EXE_SFX): test_version.o $(YLAPACK_HEADERS)
	$(CC) $(LDFLAGS) -o "$@" "$<" $(PKG_DEPLIBS) $(MATHLIB)

release: $(RELEASE_NAME)

$(RELEASE_NAME):
	@if test "x$(RELEASE_VERSION)" = "x"; then \
	  echo >&2 "set package version:  make RELEASE_VERSION=... archive"; \
	else \
          dir=`basename "$(RELEASE_NAME)" .tar.bz2`; \
	  if test "x$$dir" = "x" -o "x$$dir" = "x."; then \
	    echo >&2 "bad directory name for archive"; \
	  elif test -d "$$dir"; then \
	    echo >&2 "directory $$dir already exists"; \
	  else \
	    mkdir -p "$$dir"; \
	    cp -a $(RELEASE_FILES) "$$dir/."; \
	    rm -if "$$dir"/*~ "$$dir"/*/*~; \
	    echo "$(RELEASE_VERSION)" > "$$dir/VERSION"; \
	    tar jcf "$(RELEASE_NAME)" "$$dir"; \
	    rm -rf "$$dir"; \
	    echo "$(RELEASE_NAME) created"; \
	  fi; \
	fi;

.PHONY: clean release

# -------------------------------------------------------- end of Makefile
