dnl $Id$

dnl
dnl Michele Martone's proposal for a transition towards the GNU build system for the PSBLAS library
dnl

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
AC_INIT([MLD2P4],1.0, salvatore.filippone@uniroma2.it)

# VERSION is the file containing the PSBLAS version code
# FIXME
mld2p4_cv_version="`cat VERSION`"

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
AC_MSG_CHECKING([for PSBLAS build dir])
case $PSBLAS_DIR in 
     /*) ;; 
     *) AC_MSG_ERROR([The PSBLAS build dir must be an absolute pathname
which must be specified  --with-psblas=/path/to/psblas])
esac
if test ! -d "$PSBLAS_DIR" ; then 
   AC_MSG_ERROR([Could not find PSBLAS build dir $PSBLAS_DIR!])	
fi
AC_MSG_RESULT([$PSBLAS_DIR])
pac_cv_status_file="$PSBLAS_DIR/config.status"
if test ! -f "$pac_cv_status_file" ; then 
   AC_MSG_ERROR([Could not find configure output in $PSBLAS_DIR. 
You must first run the configure script in $PSBLAS_DIR $pac_cv_status_file.])		
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
   * ) eval "INSTALL_DIR=/usr/local/psblas";;
esac
case $libdir in 
   \/* )   eval "INSTALL_LIBDIR=$libdir";;
   * ) eval "INSTALL_LIBDIR=$INSTALL_DIR/lib";;
esac
case $includedir in 
   \/* )   eval "INSTALL_INCLUDEDIR=$includedir";;
   * ) eval "INSTALL_INCLUDEDIR=$INSTALL_DIR/include";;
esac
AC_MSG_RESULT([$INSTALL_DIR $INSTALL_INCLUDEDIR $INSTALL_LIBDIR])

###############################################################################
# Compilers detection: FC,F77,CC should be set, if found.
###############################################################################

AC_PROG_CC([xlc pgcc icc gcc ])
dnl AC_PROG_CXX

if test "X$CC" == "X" ; then
	AC_MSG_ERROR([Problem : No C compiler specified nor found!])
fi
PAC_ARG_WITH_LIBS
PAC_ARG_WITH_FLAGS(clibs,CLIBS)
PAC_ARG_WITH_FLAGS(flibs,FLIBS)
PAC_ARG_WITH_FLAGS(library-path,LIBRARYPATH)
PAC_ARG_WITH_FLAGS(include-path,INCLUDEPATH)
PAC_ARG_WITH_FLAGS(module-path,MODULE_PATH)

###############################################################################
#		BLAS  library presence checks
###############################################################################

# Note : The libmkl.a (Intel Math Kernel Library) library could be used, too.
# It is sufficient to specify it as -lmkl in the CLIBS or FLIBS or LIBS 
# and specify its path adjusting -L/path in CFLAGS.

# Right now it is a matter of user's taste when linking custom applications.
# But PSBLAS examples could take advantage of these libraries, too.


###############################################################################
#		BLACS library presence checks
###############################################################################

###############################################################################
#		METIS, SuperLU, SuperLU_Dist UMFPack libraries presence checks
###############################################################################

# Note : also SuperLUStat,superlu_gridexit,.. 
# Note : SuperLU_Dist produces superlu.a by default..
# Note : THESE WERE NOT TESTED
dnl AC_MSG_CHECKING([SuperLU_Dist])
dnl AC_CHECK_LIB(superlu,[superlu_malloc_dist],mld2p4_cv_have_superludist=yes,mld2p4_cv_have_superludist=no)
dnl AC_MSG_RESULT([$mld2p4_cv_have_superludist])
dnl AC_MSG_CHECKING([SuperLU])
dnl AC_CHECK_LIB(superlu,[superlu_malloc],mld2p4_cv_have_superlu=yes,mld2p4_cv_have_superlu=no)
dnl AC_MSG_RESULT([$mld2p4_cv_have_superlu])

PAC_CHECK_UMFPACK
if test "x$mld2p4_cv_have_umfpack" == "xyes" ; then 
   UMF_FLAGS="-DHave_UMF_ $UMF_INCLUDES"
else
   UMF_FLAGS=""
fi

PAC_CHECK_SUPERLU
if test "x$mld2p4_cv_have_superlu" == "xyes" ; then 
   SLU_FLAGS="-DHave_SLU_ $SLU_INCLUDES"
else
   SLU_FLAGS=""
fi

PAC_CHECK_SUPERLUDIST
if test "x$mld2p4_cv_have_superludist" == "xyes" ; then 
   SLUDIST_FLAGS="-DHave_SLUDist_ $SLUDIST_INCLUDES"
else
   SLUDIST_FLAGS=""
fi


# Note  : We would like to detect PSBLAS, but this is complicated by the
#         module symbols mangling rules, which are compiler specific !
#
#	  Moreover, the PSBLAS doesn't have an installer, currently.


###############################################################################
# Variable substitutions : the Make.inc.in will have these @VARIABLES@
# substituted.

AC_SUBST(PSBLAS_DIR)

AC_SUBST(INSTALL)
AC_SUBST(INSTALL_DATA)
AC_SUBST(INSTALL_DIR)
AC_SUBST(INSTALL_LIBDIR)
AC_SUBST(INSTALL_INCLUDEDIR)

AC_SUBST(SLU_FLAGS)
AC_SUBST(SLU_LIBS)
AC_SUBST(UMF_FLAGS)
AC_SUBST(UMF_LIBS)
AC_SUBST(SLUDIST_FLAGS)
AC_SUBST(SLUDIST_LIBS)


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
	SuperLU detected      : ${mld2p4_cv_have_superlu}
	SuperLU_Dist detected : ${mld2p4_cv_have_superludist}
	UMFPack detected      : ${mld2p4_cv_have_umfpack}

	If you are satisfied, run 'make' to build ${PACKAGE_NAME} and its documentation; otherwise
	type  ./configure --help=short for a complete list of configure options specific to ${PACKAGE_NAME}.
dnl	To install the program and its documentation, run 'make install' if you are root,
dnl	or run 'su -c "make install"' if you are not root.
])

###############################################################################

