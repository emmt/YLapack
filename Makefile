# where are the sources? (automatically filled in by configure script)
srcdir=.

# these values filled in by `yorick -batch make.i`
Y_MAKEDIR=
Y_EXE=
Y_EXE_PKGS=
Y_EXE_HOME=
Y_EXE_SITE=
Y_HOME_PKG=

# ----------------------------------------------------- optimization flags

# options for make command line, e.g.-   make COPT=-g TGT=exe
COPT=$(COPT_DEFAULT)
TGT=$(DEFAULT_TGT)

# ------------------------------------------------ macros for this package

PKG_NAME=ylapack
PKG_I=$(srcdir)/lapack.i

OBJS = ylapack.o

# Makefile variables PKG_CFLAGS and PKG_DEPLIBS must be defined to compile
# and link with your preferred version of LAPACK.
#
#   PKG_CFLAGS is for pre-processor directives
#
#   PKG_DEPLIBS is for libraries
#
# in PKG_CFLAGS, you may define some pre-processor macros:
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
# PKG_CFLAGS = $(MKL_DEFS)
# PKG_DEPLIBS = -L$(MKL_DIR)/lib/$(MKL_ARCH) -l$(MKL_INTERFACE) \
#    -l$(MKL_THREAD) -l$(MKL_CORE)

# To use MKL for intel64 processor with 32-bit integers:
# MKL_DIR = /opt/intel/mkl
# MKL_ARCH = intel64
# MKL_INTERFACE = mkl_gf_lp64
# MKL_THREAD = mkl_gnu_thread
# MKL_CORE = mkl_core
# MKL_DEFS = -DUSE_MKL -DINTEGER=int -I$(MKL_DIR)/include
# PKG_CFLAGS = $(MKL_DEFS)
# PKG_DEPLIBS = -L$(MKL_DIR)/lib/$(MKL_ARCH) -l$(MKL_INTERFACE) \
#    -l$(MKL_THREAD) -l$(MKL_CORE)

# To use MKL for intel64 processor with 64-bit integers:
# MKL_DIR = /opt/intel/mkl
# MKL_ARCH = intel64
# MKL_INTERFACE = mkl_gf_ilp64
# MKL_THREAD = mkl_gnu_thread
# MKL_CORE = mkl_core
# MKL_DEFS = -DUSE_MKL -DMKL_ILP64 -DINTEGER=long -I$(MKL_DIR)/include
# PKG_CFLAGS = $(MKL_DEFS)
# PKG_DEPLIBS = -L$(MKL_DIR)/lib/$(MKL_ARCH) -l$(MKL_INTERFACE) \
#    -l$(MKL_THREAD) -l$(MKL_CORE)


YLAPACK_HEADERS = \
    $(srcdir)/ylapack.h \
    $(srcdir)/ylapack_lapack.h \
    $(srcdir)/ylapack_blas.h \
    $(srcdir)/ylapack_cblas.h

YLAPACK_RULES = \
    $(srcdir)/rules/Make.atlas \
    $(srcdir)/rules/Make.default \
    $(srcdir)/rules/Make.gotoblas2 \
    $(srcdir)/rules/Make.icc \
    $(srcdir)/rules/Make.mkl \
    $(srcdir)/rules/Make.mkl_dynamic \
    $(srcdir)/rules/Make.mkl_gnu_ia32 \
    $(srcdir)/rules/Make.mkl_gnu_ilp64 \
    $(srcdir)/rules/Make.mkl_gnu_lp64 \
    $(srcdir)/rules/Make.mkl_icc_ia32 \
    $(srcdir)/rules/Make.mkl_icc_ilp64 \
    $(srcdir)/rules/Make.mkl_icc_lp64 \
    $(srcdir)/rules/Make.mkl_static \
    $(srcdir)/rules/Make.openblas

# change to give the executable a name other than yorick
PKG_EXENAME=yorick

# PKG_DEPLIBS=-Lsomedir -lsomelib   for dependencies of this package
PKG_DEPLIBS=
# set compiler (or rarely loader) flags specific to this package
PKG_CFLAGS=
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

RELEASE_FILES = \
    $(srcdir)/AUTHORS \
    $(srcdir)/LICENSE \
    $(srcdir)/Makefile \
    $(srcdir)/NEWS \
    $(srcdir)/README \
    $(srcdir)/TODO \
    $(PKG_I) $(PKG_I_EXTRA) \
    $(YLAPACK_HEADERS) \
    $(srcdir)/ylapack.c \
    $(srcdir)/lapack-test.i \
    $(YLAPACK_RULES)

RELEASE_NAME = $(PKG_NAME)-$(RELEASE_VERSION).tar.bz2

# -------------------------------- standard targets and rules (in Makepkg)

# set macros Makepkg uses in target and dependency names
# DLL_TARGETS, LIB_TARGETS, EXE_TARGETS
# are any additional targets (defined below) prerequisite to
# the plugin library, archive library, and executable, respectively
PKG_I_DEPS=$(PKG_I)
Y_DISTMAKE=distmake

ifeq (,$(strip $(Y_MAKEDIR)))
$(info *** WARNING: Y_MAKEDIR not defined, you may run 'configure' or 'yorick -batch make.i' first)
else
include $(Y_MAKEDIR)/Make.cfg
include $(Y_MAKEDIR)/Makepkg
include $(Y_MAKEDIR)/Make$(TGT)
endif

ifeq (,$(strip $(MODEL)))
else
include $(srcdir)/rules/Make.$(MODEL)
endif

# override macros Makepkg sets for rules and other macros
# Y_HOME and Y_SITE in Make.cfg may not be correct (e.g.- relocatable)
Y_HOME=$(Y_EXE_HOME)
Y_SITE=$(Y_EXE_SITE)

# reduce chance of yorick-1.5 corrupting this Makefile
MAKE_TEMPLATE = protect-against-1.5

# ------------------------------------- targets and rules for this package

# Dummy default target in case Y_MAKEDIR was not defined:
dummy-default:
	@echo >&2 "*** ERROR: Y_MAKEDIR not defined, aborting..."; false

%.o: ${srcdir}/%.c
	$(CC) -I$(srcdir) $(CPPFLAGS) $(CFLAGS) -o "$@" -c "$<"

# simple example:
#myfunc.o: myapi.h
# more complex example (also consider using PKG_CFLAGS above):
#myfunc.o: myapi.h myfunc.c
#	$(CC) $(CPPFLAGS) $(CFLAGS) -DMY_SWITCH -o $@ -c myfunc.c

FETCH_DEFS=tclsh $(srcdir)/fetch.tcl
LAPACK_LIBRARY=/usr/lib/atlas-base/atlas/liblapack.a
BLAS_LIBRARY=/usr/lib/atlas-base/atlas/libblas.a
ylapack_lapack.def: fetch.tcl
	$(FETCH_DEFS) --lapack "$(LAPACK_LIBRARY)" > "$@"
ylapack_blas.def: fetch.tcl
	$(FETCH_DEFS) --lapack "$(BLAS_LIBRARY)" > "$@"
ylapack.o: $(srcdir)/ylapack.c Makefile $(YLAPACK_HEADERS)
test_version.o: $(srcdir)/test_version.c Makefile $(YLAPACK_HEADERS)
test_version$(EXE_SFX): test_version.o $(YLAPACK_HEADERS)
	$(CC) $(LDFLAGS) -o "$@" "$<" $(PKG_DEPLIBS) $(MATHLIB)

release: $(RELEASE_NAME)

$(RELEASE_NAME):
	@if test "x$(RELEASE_VERSION)" = "x"; then \
	  echo >&2 "set package version:  make RELEASE_VERSION=... release"; \
	else \
          dir=`basename "$(RELEASE_NAME)" .tar.bz2`; \
	  if test "x$$dir" = "x" -o "x$$dir" = "x."; then \
	    echo >&2 "bad directory name for archive"; \
	  elif test -d "$$dir"; then \
	    echo >&2 "directory $$dir already exists"; \
	  else \
	    mkdir -p "$$dir" "$$dir/rules"; \
	    for src in $(RELEASE_FILES); do \
	      base=`basename "$$src"`; \
	      if echo "x$$src" | grep -q '/rules/[^/]*$$'; then \
	        dst="$$dir/rules/$$base"; \
	      else \
	        dst="$$dir/$$base"; \
	      fi; \
	      if test "$$base" = "Makefile"; then \
	        sed <"$$src" >"$$dst" -e 's/^\( *Y_\(MAKEDIR\|EXE\(\|_PKGS\|_HOME\|_SITE\)\|HOME_PKG\) *=\).*/\1/'; \
	        touch -r "$$src" "$$dst"; \
	      else \
	        cp -p "$$src" "$$dst"; \
	      fi; \
	    done; \
	    rm -f "$$dir"/*~ "$$dir"/*/*~; \
	    echo "$(RELEASE_VERSION)" > "$$dir/VERSION"; \
	    tar jcf "$(RELEASE_NAME)" "$$dir"; \
	    rm -rf "$$dir"; \
	    echo "$(RELEASE_NAME) created"; \
	  fi; \
	fi;

.PHONY: clean release

# -------------------------------------------------------- end of Makefile
