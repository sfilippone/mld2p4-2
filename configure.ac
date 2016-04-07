dnl $Id$

dnl Process this file with autoconf to produce a configure script.
dnl
dnl usage : aclocal -I config/ && autoconf && ./configure && make 

dnl then  : VAR=VAL ./configure

dnl In some configurations (AIX) the next line is needed:
dnl MPIFC=mpxlf95 ./configure

dnl then  : ./configure VAR=VAL 
dnl then  : ./configure --help=short
dnl then  : ./configure --help

dnl the PSBLAS modules get this task difficult to accomplish!
dnl SEE  : --module-path --include-path

dnl NOTE : There is no cross compilation support.

dnl NOTE : missing ifort and kl* library handling..
dnl NOTE : odd configurations like ifc + gcc still await in the mist of the unknown


###############################################################################
###############################################################################
# 
#  This script is used by the PSBLAS to determine the compilers, linkers, and
# libraries to build its libraries executable code.
#  Its behaviour is driven on the compiler it finds or it is dictated to work
# with.
#
###############################################################################
###############################################################################

# NOTE: the literal for version (the second argument to AC_INIT should be a literal!)
AC_INIT([MLD2P4],2.0, bugreport@mld2p4.it)

# VERSION is the file containing the PSBLAS version code
# FIXME
mld2p4_cv_version="2.0"

# A sample source file
AC_CONFIG_SRCDIR([mlprec/mld_prec_type.f90])
# Our custom M4 macros are in the 'config' directory
AC_CONFIG_MACRO_DIR([config])
AC_MSG_NOTICE([
--------------------------------------------------------------------------------
	Welcome to the $PACKAGE_NAME $mld2p4_cv_version configure Script.

	This creates Make.inc, but if you read carefully the
	documentation, you can make your own by hand for your needs.

        ./configure  --with-psblas=/path/to/psblas
	See ./configure --help=short fore more info.
--------------------------------------------------------------------------------
		])

###############################################################################
# 			FLAGS and LIBS user customization
###############################################################################

dnl NOTE : no spaces before the comma, and no brackets before the second argument!
PAC_ARG_WITH_PSBLAS
PSBLAS_DIR="$pac_cv_psblas_dir";
PSBLAS_INCDIR="$pac_cv_psblas_incdir";
PSBLAS_LIBDIR="$pac_cv_psblas_libdir";
AC_MSG_CHECKING([for PSBLAS install dir])
case $PSBLAS_DIR in 
     /*) ;; 
     *) AC_MSG_ERROR([The PSBLAS installation dir must be an absolute pathname
 specified with --with-psblas=/path/to/psblas])
esac
if test ! -d "$PSBLAS_DIR" ; then 
   AC_MSG_ERROR([Could not find PSBLAS build dir $PSBLAS_DIR!])	
fi
AC_MSG_RESULT([$PSBLAS_DIR])
pac_cv_status_file="$PSBLAS_INCDIR/Make.inc.psblas"
if test ! -f "$pac_cv_status_file" ; then 
   AC_MSG_ERROR([Could not find an installation  in $PSBLAS_DIR.])		
fi
dnl. $pac_cv_status_file
dnl AC_MSG_NOTICE([Loaded $pac_cv_status_file $FC $MPIFC $BLACS_LIBS])
dnl AM_INIT_AUTOMAKE
dnl Specify required version of autoconf.
AC_PREREQ(2.59)
#
# Installation. 
#
#
AC_PROG_INSTALL

AC_MSG_CHECKING([where to install])
case $prefix in  
   \/* )   eval "INSTALL_DIR=$prefix";;
   * ) eval "INSTALL_DIR=/usr/local/mld2p4";;
esac
case $libdir in 
   \/* )   eval "INSTALL_LIBDIR=$libdir";;
   * ) eval "INSTALL_LIBDIR=$INSTALL_DIR/lib";;
esac
case $includedir in 
   \/* )   eval "INSTALL_INCLUDEDIR=$includedir";;
   * ) eval "INSTALL_INCLUDEDIR=$INSTALL_DIR/include";;
esac
case $docsdir in 
   \/* )   eval "INSTALL_DOCSDIR=$docsdir";;
   * ) eval "INSTALL_DOCSDIR=$INSTALL_DIR/docs";;
esac
AC_MSG_RESULT([$INSTALL_DIR $INSTALL_INCLUDEDIR $INSTALL_LIBDIR $INSTALL_DOCSDIR])

###############################################################################
# Compilers detection: FC,F77,CC should be set, if found.
###############################################################################

AC_PROG_FC([ftn xlf95 xlf90 xlf pgf95 pgf90 ifort ifc nagfor  gfortran])
AC_PROG_CC([xlc pgcc icc gcc cc])
dnl AC_PROG_CXX

if test "X$CC" == "X" ; then
	AC_MSG_ERROR([Problem : No C compiler specified nor found!])
fi
if eval "$FC -qversion 2>&1 | grep XL 2>/dev/null" ; then
	# Some configurations of the XLF want "-WF," prepended to -D.. flags.
	# TODO : discover the exact conditions when the usage of -WF is needed.
	mld_cv_define_prepend="-WF,"
        if eval "$MPIFC -qversion 2>&1 | grep -e\"Version: 10\.\" 2>/dev/null"; then
        	FDEFINES="$mld_cv_define_prepend-DXLF_10 $FDEFINES"		  
        fi

 # Note : there coule be problems with old xlf compiler versions ( <10.1 )
 # since (as far as it is known to us) -WF, is not used in earlier versions.
 # More problems could be undocumented yet.
fi 
PAC_ARG_WITH_EXTRA_LIBS
dnl PAC_ARG_WITH_FLAGS(clibs,CLIBS)
dnl PAC_ARG_WITH_FLAGS(flibs,FLIBS)
dnl PAC_ARG_WITH_FLAGS(library-path,LIBRARYPATH)
dnl PAC_ARG_WITH_FLAGS(include-path,INCLUDEPATH)
dnl PAC_ARG_WITH_FLAGS(module-path,MODULE_PATH)

dnl CPPFLAGS="$INCLUDEPATH $CPPFLAGS"
###############################################################################
#		BLAS  library presence checks
###############################################################################

# Note : The libmkl.a (Intel Math Kernel Library) library could be used, too.
# It is sufficient to specify it as -lmkl in the CLIBS or FLIBS or LIBS 
# and specify its path adjusting -L/path in CFLAGS.

# Right now it is a matter of user's taste when linking custom applications.
# But PSBLAS examples could take advantage of these libraries, too.


###############################################################################
#	
###############################################################################
PAC_ARG_SERIAL_MPI

#Note : we miss the name of the Intel C compiler
if test x"$pac_cv_serial_mpi" == x"yes" ; then
   FAKEMPI="fakempi.o";
   MPIFC="$FC";
   MPIF77="$F77";
   MPICC="$CC";
else 
AC_LANG([C])
if test "X$MPICC" = "X" ; then
    # This is our MPICC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPICC],[mpxlc mpcc mpicc cc])
fi
ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for C]])])
AC_LANG(Fortran 77)
if test "X$MPIF77" = "X" ; then
    # This is our MPIFC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPIF77],[mpxlf mpf77 mpif77 ftn])
fi
ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for Fortran 77]])])

AC_LANG([Fortran])

if test "X$MPIFC" = "X" ; then
    # This is our MPIFC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPIFC],[mpxlf2003_r mpxlf2003 mpxlf95_r mpxlf90 mpf95 mpf90  mpif95 mpif90  ])
fi

ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for Fortran]])])
fi

# We leave a default language for the next checks.
dnl AC_LANG([Fortran 77])
AC_LANG([C])

dnl Now on, MPIFC should be set, as MPIF77 and MPICC

###############################################################################
# Sanity checks, although redundant (useful when debugging this configure.ac)!
###############################################################################

###############################################################################
#		PSBLAS library presence checks
###############################################################################
AX_F90_MODULE_FLAG
FMFLAG="${ax_cv_f90_modflag%%[ ]*}"
# Last resort attempt
if test "x$FMFLAG" == "xnot" ; then 
   AC_MSG_NOTICE([Fortran inclusion flag detection failed, trying with -I])  
   MODEXT=".mod"
   FMFLAG="-I"
   FIFLAG="-I"
fi


PAC_FORTRAN_HAVE_PSBLAS([AC_MSG_RESULT([yes.])],
	[AC_MSG_ERROR([no.  Could not find working version of PSBLAS.])])

PAC_FORTRAN_PSBLAS_VERSION()
if test "x$pac_cv_psblas_major" == "xunknown"; then
  AC_MSG_ERROR([PSBLAS version major "$pac_cv_psblas_major".])
fi
if test "x$pac_cv_psblas_minor" == "xunknown"; then
   AC_MSG_ERROR([PSBLAS version minor "$pac_cv_psblas_minor".])
fi
if test "x$pac_cv_psblas_patchlevel" == "xunknown"; then
   AC_MSG_ERROR([PSBLAS patchlevel "$pac_cv_psblas_patchlevel".])
fi
if (( $pac_cv_psblas_major < 3 )) ||
   (  (( $pac_cv_psblas_major == 3 )) && (( $pac_cv_psblas_minor < 4 ))) ; then
  AC_MSG_ERROR([I need at least PSBLAS version 3.4.]) 
else
  AC_MSG_NOTICE([Am configuring with PSBLAS version $pac_cv_psblas_major.$pac_cv_psblas_minor.$pac_cv_psblas_patchlevel.])
fi
###############################################################################
#		Auxiliary packages
###############################################################################


PAC_CHECK_MUMPS
if test "x$mld2p4_cv_have_mumps" == "xyes" ; then 
   MUMPS_FLAGS="-DHave_MUMPS_ $MUMPS_INCLUDES"
   FDEFINES="$mld_cv_define_prepend-DHAVE_MUMPS_ $mld2p4_cv_mumpsincdir $FDEFINES"
else
   MUMPS_FLAGS=""
fi

PAC_CHECK_UMFPACK
if test "x$mld2p4_cv_have_umfpack" == "xyes" ; then 
   UMF_FLAGS="-DHave_UMF_ $UMF_INCLUDES"
   FDEFINES="$mld_cv_define_prepend-DHAVE_UMF_ $FDEFINES"
else
   UMF_FLAGS=""
fi

PAC_CHECK_SUPERLU
if test "x$mld2p4_cv_have_superlu" == "xyes" ; then 
   SLU_FLAGS="-DHave_SLU_ -DSLU_VERSION_$pac_slu_version $SLU_INCLUDES"
   FDEFINES="$mld_cv_define_prepend-DHAVE_SLU_ $FDEFINES"
else
   SLU_FLAGS=""
fi

PAC_CHECK_SUPERLUDIST
if test "x$mld2p4_cv_have_superludist" == "xyes" ; then 
   SLUDIST_FLAGS=""	       
   SLUDIST_FLAGS="-DHave_SLUDist_ -DSLUD_VERSION_$pac_sludist_version $SLUDIST_INCLUDES"
   FDEFINES="$mld_cv_define_prepend-DHAVE_SLUDIST_ $FDEFINES"
else
   SLUDIST_FLAGS=""
fi


###############################################################################
# Variable substitutions : the Make.inc.in will have these @VARIABLES@
# substituted.

AC_SUBST(PSBLAS_DIR)
AC_SUBST(PSBLAS_INCDIR)

AC_SUBST(INSTALL)
AC_SUBST(INSTALL_DATA)
AC_SUBST(INSTALL_DIR)
AC_SUBST(INSTALL_LIBDIR)
AC_SUBST(INSTALL_INCLUDEDIR)
AC_SUBST(INSTALL_DOCSDIR)

AC_SUBST(EXTRA_LIBS)
AC_SUBST(MUMPS_FLAGS)
AC_SUBST(MUMPS_LIBS)
AC_SUBST(SLU_FLAGS)
AC_SUBST(SLU_LIBS)
AC_SUBST(UMF_FLAGS)
AC_SUBST(UMF_LIBS)
AC_SUBST(SLUDIST_FLAGS)
AC_SUBST(SLUDIST_LIBS)
AC_SUBST(FDEFINES)

###############################################################################
# the following files will be created by Automake

AC_CONFIG_FILES([Make.inc])
AC_OUTPUT()
#AC_OUTPUT(Make.inc Makefile)
###############################################################################

dnl Please note that brackets around variable identifiers are absolutely needed for compatibility..
AC_MSG_NOTICE([
	${PACKAGE_NAME} ${mld2p4_cv_version} has been configured as follows:

	PSBLAS library        : ${PSBLAS_DIR}
	MUMPS detected        : ${mld2p4_cv_have_mumps}
	SuperLU detected      : ${mld2p4_cv_have_superlu}
	SuperLU_Dist detected : ${mld2p4_cv_have_superludist}
	UMFPack detected      : ${mld2p4_cv_have_umfpack}

	If you are satisfied, run 'make' to build ${PACKAGE_NAME} and its documentation; otherwise
	type  ./configure --help=short for a complete list of configure options specific to ${PACKAGE_NAME}.
dnl	To install the program and its documentation, run 'make install' if you are root,
dnl	or run 'su -c "make install"' if you are not root.
])

###############################################################################

