dnl
dnl $Id$
dnl
dnl 20080206
dnl M4 macros for the PSBLAS library and useful for packages using PSBLAS.
dnl

dnl @synopsis PAC_CHECK_LIBS
dnl
dnl Tries to detect the presence of a specific function among various libraries, using AC_CHECK_LIB
dnl repeatedly on the specified libraries.
dnl 
dnl Example use:
dnl
dnl PAC_CHECK_LIBS([atlas blas],
dnl		[dgemm],
dnl		[have_dgemm=yes],
dnl		[have_dgemm=no])
dnl 
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
dnl 20080211 modified slighty from original.
AC_DEFUN([PAC_CHECK_LIBS],
[
 pac_check_libs_ok=no
 [for pac_check_libs_f in $2 
 do ]
 [for pac_check_libs_l in $1 
 do ]
    if test x"$pac_check_libs_ok" == xno ; then
     AC_CHECK_LIB([$pac_check_libs_l],[$pac_check_libs_f], [pac_check_libs_ok=yes; pac_check_libs_LIBS="-l$pac_check_libs_l"],[],[$5])
    fi
  done
  done
 # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
 [ if test x"$pac_check_libs_ok" = xyes ; then
	$3
 else
        pac_check_libs_ok=no
        $4
 fi
 ]
])dnl 


dnl @synopsis PAC_FORTRAN_HAVE_MOVE_ALLOC( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program with move_alloc (a Fortran 2003 function).
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN([PAC_FORTRAN_HAVE_MOVE_ALLOC],
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran MOVE_ALLOC intrinsic])
 AC_LANG_PUSH([Fortran])
 ac_ext='f90';
 AC_COMPILE_IFELSE([ program test_move_alloc
		       integer, allocatable :: a(:), b(:)
		       allocate(a(3))
		       call move_alloc(a, b)
		       print *, allocated(a), allocated(b)
		       print *, b
		     end program test_move_alloc],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])


dnl @synopsis PAC_CHECK_HAVE_GFORTRAN( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will check if MPIFC is $FC.
dnl The check will proceed by compiling a small Fortran program
dnl containing the __GNUC__ macro, which should be defined in the
dnl gfortran compiled programs.
dnl
dnl On pass, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_HAVE_GFORTRAN,
[AC_MSG_CHECKING([for GNU Fortran])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program main
#ifdef __GNUC__ 
              print *, "GCC!"
#else
        this program will fail
#endif
           end],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])


dnl @synopsis PAC_HAVE_MODERN_GFORTRAN( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will check if the GNU fortran version is suitable for PSBLAS.
dnl If yes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl Note : Will use MPIFC; if unset, will use '$FC'.
dnl 
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN(PAC_HAVE_MODERN_GFORTRAN,
 [AC_MSG_CHECKING([for recent GNU Fortran])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program main
#if ( __GNUC__ >= 4 && __GNUC_MINOR__ >= 6 ) || ( __GNUC__ > 4 )
              print *, "ok"
#else
        this program will fail
#endif
           end],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_CHECK_HAVE_MPI_MOD( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will determine if the fortran compiler MPIFC needs to include mpi.h or needs
dnl to use the mpi module.
dnl
dnl If yes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl 
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl Modified Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_CHECK_HAVE_MPI_MOD,
 [AC_MSG_CHECKING([for Fortran MPI mod])
  AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program test
             use mpi
           end program test],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_ARG_WITH_FLAGS(lcase_name, UCASE_NAME)
dnl
dnl Test for --with-lcase_name="compiler/loader flags".  if defined, prepend 
dnl flags to standard UCASE_NAME definition.
dnl
dnl Use this macro to facilitate additional special flags that should be
dnl passed on to the preprocessor/compilers/loader.
dnl
dnl NOTE : Renamed after TAC_ARG_WITH_FLAGS as in the Trilinos-8.0.4 package.
dnl 
dnl NOTE : This macro works in a way the user should invoke
dnl         --with-flags=...
dnl	   only once, otherwise the first one will take effect.
dnl
dnl Example use:
dnl 
dnl PAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
dnl 
dnl tests for --with-cxxflags and pre-pends to CXXFLAGS
dnl 
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl @notes  Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_WITH_FLAGS],
[
AC_MSG_CHECKING([whether additional [$2] flags should be added (should be invoked only once)])
dnl AC_MSG_CHECKING([whether additional [$2] flags should be added])
AC_ARG_WITH($1,
AC_HELP_STRING([--with-$1], 
[additional [$2] flags to be added: will prepend to [$2]]),
[
$2="${withval} ${$2}"
AC_MSG_RESULT([$2 = ${$2}])
],
AC_MSG_RESULT(no)
)
])


dnl @synopsis PAC_ARG_WITH_LIBS
dnl
dnl Test for --with-libs="name(s)".
dnl 
dnl Prepends the specified name(s) to the list of libraries to link 
dnl with.  
dnl
dnl note: Renamed after PAC_ARG_WITH_LIBS as in the Trilinos package.
dnl
dnl Example use:
dnl
dnl PAC_ARG_WITH_LIBS
dnl 
dnl tests for --with-libs and pre-pends to LIBS
dnl
dnl @author Jim Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([PAC_ARG_WITH_LIBS],
[
AC_MSG_CHECKING([whether additional libraries are needed])
AC_ARG_WITH(libs,
AC_HELP_STRING([--with-libs], 
[List additional link flags  here.  For example, --with-libs=-lspecial_system_lib
or --with-libs=-L/path/to/libs]),
[
LIBS="${withval} ${LIBS}"
AC_MSG_RESULT([LIBS = ${LIBS}])
],
AC_MSG_RESULT(no)
)
]
)


dnl @synopsis PAC_ARG_WITH_EXTRA_LIBS
dnl
dnl Test for --with-extra-libs="name(s)".
dnl 
dnl Appends the specified name(s) to the list of libraries to link 
dnl with.  
dnl
dnl note: Renamed after PAC_ARG_WITH_EXTRA_LIBS as in the Trilinos package.
dnl
dnl Example use:
dnl
dnl PAC_ARG_WITH_EXTRA_LIBS
dnl 
dnl tests for --with-extra-libs and pre-pends to LIBS
dnl
dnl @author Jim Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([PAC_ARG_WITH_EXTRA_LIBS],
[
AC_MSG_CHECKING([whether additional libraries are needed])
AC_ARG_WITH(extra-libs,
AC_HELP_STRING([--with-extra-libs], 
[List additional link flags  here.  For example, --with-extra-libs=-lspecial_system_lib
or --with-extra-libs=-L/path/to/libs]),
[
EXTRA_LIBS="${withval}"
AC_MSG_RESULT([EXTRA_LIBS = ${EXTRA_LIBS}])
],
AC_MSG_RESULT(no)
)
]
)

dnl @synopsis PAC_ARG_WITH_PSBLAS
dnl
dnl Test for --with-psblas="pathname".
dnl 
dnl Defines the path to PSBLAS build dir.
dnl
dnl note: Renamed after PAC_ARG_WITH_LIBS as in the Trilinos package.
dnl
dnl Example use:
dnl
dnl PAC_ARG_WIT_PSBLAS
dnl 
dnl tests for --with-psblas and pre-pends to PSBLAS_PATH
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_WITH_PSBLAS],
[
AC_ARG_WITH(psblas,
AC_HELP_STRING([--with-psblas=DIR], [The install directory for PSBLAS, for example,
 --with-psblas=/opt/packages/psblas-3.3]),
[pac_cv_psblas_dir=$withval],
[pac_cv_psblas_dir=''])
AC_ARG_WITH(psblas-incdir, AC_HELP_STRING([--with-psblas-incdir=DIR], [Specify the directory for PSBLAS includes.]),
        [pac_cv_psblas_incdir=$withval],
        [pac_cv_psblas_incdir=''])
AC_ARG_WITH(psblas-libdir, AC_HELP_STRING([--with-psblas-libdir=DIR], [Specify the directory for PSBLAS library.]),
        [pac_cv_psblas_libdir=$withval],
        [pac_cv_psblas_libdir=''])
if test x"$pac_cv_psblas_incdir" == "x" ; then
   pac_cv_psblas_incdir="$pac_cv_psblas_dir/include";
fi
if test x"$pac_cv_psblas_libdir" == "x" ; then
   pac_cv_psblas_libdir="$pac_cv_psblas_dir/lib";
fi

])

dnl @synopsis PAC_FORTRAN_HAVE_PSBLAS( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program using the PSBLAS library
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_HAVE_PSBLAS,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([for working installation of PSBLAS])
 AC_LANG_PUSH([Fortran])
 ac_objext='o'
 ac_ext='f90'
 ac_fc="${MPIFC-$FC}";
 save_FCFLAGS="$FCFLAGS";
 FCFLAGS=" $FMFLAG$PSBLAS_DIR/include $save_FCFLAGS"
 AC_COMPILE_IFELSE([
		    program test
		    use psb_base_mod
		    end program test],
		   [ifelse([$1], , :, [ $1   ])],
		   [  ifelse([$2], , , [ $2 ])])
AC_LANG_POP([Fortran])
rm -f conftest*])


dnl @synopsis PAC_FORTRAN_PSBLAS_VERSION( )
dnl
dnl Will try to compile, link and run  a program using the PSBLAS library. \
dnl  Checks for version major,  minor and patchlevel
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_PSBLAS_VERSION,
[AC_MSG_CHECKING([for version of PSBLAS])
AC_LANG_PUSH([Fortran])	  
ac_exeext=''
ac_objext='o'
ac_ext='f90'
save_FCFLAGS=$FCFLAGS;
FCFLAGS=" $FMFLAG$PSBLAS_DIR/include $save_FCFLAGS"
save_LDFLAGS=$LDFLAGS;
LDFLAGS=" -L$PSBLAS_DIR/lib -lpsb_base $save_LDFLAGS"

dnl ac_compile='${MPIFC-$FC} -c -o conftest${ac_objext} $FMFLAG$PSBLAS_DIR/include $FMFLAG$PSBLAS_DIR/lib conftest.$ac_ext  1>&5'
dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FCFLAGS $LDFLAGS conftest.$ac_ext $FMFLAG$PSBLAS_DIR/include -L$PSBLAS_DIR/lib -lpsb_base $LIBS 1>&5'
dnl Warning : square brackets are EVIL!

AC_LINK_IFELSE([
		program test
		use psb_base_mod
		print *,psb_version_major_
		end program test],
	       [pac_cv_psblas_major=`./conftest${ac_exeext} | sed 's/^ *//'`],
	       [pac_cv_psblas_major="unknown"])
  
AC_LINK_IFELSE([
		program test
		use psb_base_mod
		print *,psb_version_minor_
		end program test],
	       [pac_cv_psblas_minor=`./conftest${ac_exeext} | sed 's/^ *//'`],
	       [pac_cv_psblas_minor="unknown"])
  
AC_LINK_IFELSE([
		program test
		use psb_base_mod
		print *,psb_patchlevel_
		end program test],
	       [pac_cv_psblas_patchlevel=`./conftest${ac_exeext} | sed 's/^ *//'`],
	       [pac_cv_psblas_patchlevel="unknown"])

AC_MSG_RESULT([Done])
AC_LANG_POP([Fortran])])

dnl @synopsis PAC_FORTRAN_TEST_TR15581( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the TR15581 Fortran extension support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_TR15581,
[AC_MSG_CHECKING([support for Fortran allocatables TR15581])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
module conftest
  type outer
    integer,  allocatable :: v(:)
  end type outer

  interface foo
    module procedure foov, food
  end interface
contains

  subroutine foov(a,b)

    implicit none
    integer, allocatable, intent(inout) :: a(:)
    integer, allocatable, intent(out) :: b(:)


    allocate(b(size(a)))

  end subroutine foov
  subroutine food(a,b)

    implicit none
    type(outer), intent(inout) :: a
    type(outer), intent(out) :: b


    allocate(b%v(size(a%v)))

  end subroutine food

end module conftest



program testtr15581
  use conftest
  type(outer) :: da, db
  integer, allocatable :: a(:), b(:)

  allocate(a(10),da%v(10))
  a = (/ (i,i=1,10) /)
  da%v = (/ (i,i=1,10) /)
  call foo(a,b)
  call foo(da,db)
  write(*,*) b
  write(*,*) db%v

end program testtr15581],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_TEST_VOLATILE( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the VOLATILE Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_VOLATILE,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran VOLATILE])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program conftest
  integer, volatile :: i, j
end program conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AC_FD_CC
		     cat conftest.$ac_ext >&AC_FD_CC
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_MAKE_IS_GNUMAKE
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
define(PAC_MAKE_IS_GNUMAKE,[
AC_MSG_CHECKING(for gnumake)
MAKE=${MAKE:-make}

if $MAKE --version 2>&1 | grep -e"GNU Make" >/dev/null; then 
    AC_MSG_RESULT(yes)
    psblas_make_gnumake='yes'
else
    AC_MSG_RESULT(no)
    psblas_make_gnumake='no'
fi
])dnl

dnl @synopsis PAC_CHECK_UMFPACK
dnl
dnl Will try to find the UMFPACK library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_UMFPACK,
[AC_ARG_WITH(umfpack, AC_HELP_STRING([--with-umfpack=LIBNAME], [Specify the library name for UMFPACK and its support libraries. 
Default: "-lumfpack -lamd"]),
        [mld2p4_cv_umfpack=$withval],
        [mld2p4_cv_umfpack='-lumfpack -lamd'])
AC_ARG_WITH(umfpackdir, AC_HELP_STRING([--with-umfpackdir=DIR], [Specify the directory for UMFPACK library and includes.]),
        [mld2p4_cv_umfpackdir=$withval],
        [mld2p4_cv_umfpackdir=''])
AC_ARG_WITH(umfpackincdir, AC_HELP_STRING([--with-umfpackincdir=DIR], [Specify the directory for UMFPACK includes.]),
        [mld2p4_cv_umfpackincdir=$withval],
        [mld2p4_cv_umfpackincdir=''])
AC_ARG_WITH(umfpacklibdir, AC_HELP_STRING([--with-umfpacklibdir=DIR], [Specify the directory for UMFPACK library.]),
        [mld2p4_cv_umfpacklibdir=$withval],
        [mld2p4_cv_umfpacklibdir=''])

AC_LANG_PUSH([C])
save_LIBS="$LIBS"
save_CPPFLAGS="$CPPFLAGS"
if test "x$mld2p4_cv_umfpackincdir" != "x"; then 
 AC_MSG_NOTICE([umfp include dir $mld2p4_cv_umfpackincdir])
 UMF_INCLUDES="-I$mld2p4_cv_umfpackincdir"
elif test "x$mld2p4_cv_umfpackdir" != "x"; then 
 AC_MSG_NOTICE([umfp dir $mld2p4_cv_umfpackdir])
 UMF_INCLUDES="-I$mld2p4_cv_umfpackdir" 
fi
CPPFLAGS="$UMF_INCLUDES $CPPFLAGS"
AC_CHECK_HEADER([umfpack.h],
 [pac_umf_header_ok=yes],
 [pac_umf_header_ok=no; UMF_INCLUDES=""])
if test "x$mld2p4_cv_umfpacklibdir" != "x"; then 
   LIBS="-L$mld2p4_cv_umfpacklibdir $LIBS"
   UMF_LIBDIR="-L$mld2p4_cv_umfpacklibdir"
elif test "x$mld2p4_cv_umfpackdir" != "x"; then 
   LIBS="-L$mld2p4_cv_umfpackdir $LIBS"
   UMF_LIBDIR="-L$mld2p4_cv_umfpackdir"
fi
if test "x$pac_umf_header_ok" == "xno" ; then
  unset ac_cv_header_umfpack_h
  UMF_INCLUDES="-I$mld2p4_cv_umfpackdir" 
  CPPFLAGS="$UMF_INCLUDES $SAVE_CPPFLAGS"
  AC_CHECK_HEADER([umfpack.h],
  [pac_umf_header_ok=yes],
  [pac_umf_header_ok=no; UMF_INCLUDES=""])
fi
if test "x$pac_umf_header_ok" == "xno" ; then
dnl Maybe Include or include subdirs? 
  unset ac_cv_header_umfpack_h
  UMF_INCLUDES="-I$mld2p4_cv_umfpackdir/include -I$mld2p4_cv_umfpackdir/Include "
  CPPFLAGS="$UMF_INCLUDES $SAVE_CPPFLAGS"

  AC_MSG_CHECKING([for umfpack_di_symbolic in $UMF_INCLUDES])
  AC_CHECK_HEADER([umfpack.h],
    [pac_umf_header_ok=yes],
    [pac_umf_header_ok=no; UMF_INCLUDES=""])
fi
if test "x$pac_umf_header_ok" == "xno" ; then
dnl Maybe new structure with UMFPACK UFconfig AMD? 
   unset ac_cv_header_umfpack_h
   UMF_INCLUDES="-I$mld2p4_cv_umfpackdir/UFconfig -I$mld2p4_cv_umfpackdir/UMFPACK/Include -I$mld2p4_cv_umfpackdir/AMD/Include"
   CPPFLAGS="$UMF_INCLUDES $SAVE_CPPFLAGS"
   AC_CHECK_HEADER([umfpack.h],
     [pac_umf_header_ok=yes],
     [pac_umf_header_ok=no; UMF_INCLUDES=""])
fi


if test "x$pac_umf_header_ok" == "xyes" ; then 
      UMF_LIBS="$mld2p4_cv_umfpack $UMF_LIBDIR"
      LIBS="$UMF_LIBS -lm $LIBS";
      AC_MSG_CHECKING([for umfpack_di_symbolic in $UMF_LIBS])
      AC_TRY_LINK_FUNC(umfpack_di_symbolic, 
       [mld2p4_cv_have_umfpack=yes;pac_umf_lib_ok=yes; ],
       [mld2p4_cv_have_umfpack=no;pac_umf_lib_ok=no; UMF_LIBS=""])
      AC_MSG_RESULT($pac_umf_lib_ok)
     if test "x$pac_umf_lib_ok" == "xno" ; then 
        dnl Maybe Lib or lib? 
        UMF_LIBDIR="-L$mld2p4_cv_umfpackdir/Lib -L$mld2p4_cv_umfpackdir/lib"
        UMF_LIBS="$mld2p4_cv_umfpack $UMF_LIBDIR -lm $SAVE_LIBS"
        LIBS="$UMF_LIBS"
        
      AC_MSG_CHECKING([for umfpack_di_symbolic in $UMF_LIBS])
      AC_TRY_LINK_FUNC(umfpack_di_symbolic, 
       [mld2p4_cv_have_umfpack=yes;pac_umf_lib_ok=yes; ],
       [mld2p4_cv_have_umfpack=no;pac_umf_lib_ok=no; UMF_LIBS=""])
      AC_MSG_RESULT($pac_umf_lib_ok)
     fi
     if test "x$pac_umf_lib_ok" == "xno" ; then 
        dnl Maybe UMFPACK/Lib? 
        UMF_LIBDIR="-L$mld2p4_cv_umfpackdir/AMD/Lib -L$mld2p4_cv_umfpackdir/UMFPACK/Lib"
        UMF_LIBS="$mld2p4_cv_umfpack $UMF_LIBDIR -lm $SAVE_LIBS"
             LIBS="$UMF_LIBS"
      AC_MSG_CHECKING([for umfpack_di_symbolic in $UMF_LIBS])
      AC_TRY_LINK_FUNC(umfpack_di_symbolic, 
       [mld2p4_cv_have_umfpack=yes;pac_umf_lib_ok=yes; ],
       [mld2p4_cv_have_umfpack=no;pac_umf_lib_ok=no; UMF_LIBS=""])
      AC_MSG_RESULT($pac_umf_lib_ok)
     fi
fi
LIBS="$SAVE_LIBS";
CPPFLAGS="$SAVE_CPPFLAGS";
AC_LANG_POP([C])
])dnl 

dnl @synopsis PAC_CHECK_SUPERLU
dnl
dnl Will try to find the SUPERLU library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_SUPERLU,
[AC_ARG_WITH(superlu, AC_HELP_STRING([--with-superlu=LIBNAME], [Specify the library name for SUPERLU library.
Default: "-lsuperlu"]),
        [mld2p4_cv_superlu=$withval],
        [mld2p4_cv_superlu='-lsuperlu'])
AC_ARG_WITH(superludir, AC_HELP_STRING([--with-superludir=DIR], [Specify the directory for SUPERLU library and includes.]),
        [mld2p4_cv_superludir=$withval],
        [mld2p4_cv_superludir=''])
AC_ARG_WITH(superluincdir, AC_HELP_STRING([--with-superluincdir=DIR], [Specify the directory for SUPERLU includes.]),
        [mld2p4_cv_superluincdir=$withval],
        [mld2p4_cv_superluincdir=''])
AC_ARG_WITH(superlulibdir, AC_HELP_STRING([--with-superlulibdir=DIR], [Specify the directory for SUPERLU library.]),
        [mld2p4_cv_superlulibdir=$withval],
        [mld2p4_cv_superlulibdir=''])
AC_LANG_PUSH([C])
save_LIBS="$LIBS"
save_CPPFLAGS="$CPPFLAGS"
if test "x$mld2p4_cv_superluincdir" != "x"; then 
 AC_MSG_NOTICE([slu include dir $mld2p4_cv_superluincdir])
 SLU_INCLUDES="-I$mld2p4_cv_superluincdir"
elif test "x$mld2p4_cv_superludir" != "x"; then 
 AC_MSG_NOTICE([slu dir $mld2p4_cv_superludir])
 SLU_INCLUDES="-I$mld2p4_cv_superludir"
fi
if test "x$mld2p4_cv_superlulibdir" != "x"; then 
   SLU_LIBS="-L$mld2p4_cv_superlulibdir"
elif test "x$mld2p4_cv_superludir" != "x"; then 
   SLU_LIBS="-L$mld2p4_cv_superludir"
fi

LIBS="$SLU_LIBS $LIBS"
CPPFLAGS="$SLU_INCLUDES $save_CPPFLAGS"
AC_CHECK_HEADERS([slu_ddefs.h],
		[pac_slu_header_ok=yes],
		[pac_slu_header_ok=no; SLU_INCLUDES=""])
if test "x$pac_slu_header_ok" == "xno" ; then 
dnl Maybe Include or include subdirs? 
  unset ac_cv_header_slu_ddefs_h
  SLU_INCLUDES="-I$mld2p4_cv_superludir/include -I$mld2p4_cv_superludir/Include "
  CPPFLAGS="$SLU_INCLUDES $save_CPPFLAGS"

  AC_CHECK_HEADERS([slu_ddefs.h],
		 [pac_slu_header_ok=yes],
		 [pac_slu_header_ok=no; SLU_INCLUDES=""])
fi

if test "x$pac_slu_header_ok" == "xyes" ; then 
 SLU_LIBS="$mld2p4_cv_superlu $SLU_LIBS"
 LIBS="$SLU_LIBS -lm $save_LIBS";
 AC_MSG_CHECKING([for superlu_malloc in $SLU_LIBS])
 AC_TRY_LINK_FUNC(superlu_malloc, 
		  [mld2p4_cv_have_superlu=yes;pac_slu_lib_ok=yes;],
		  [mld2p4_cv_have_superlu=no;pac_slu_lib_ok=no; SLU_LIBS=""; ])
 if test "x$pac_slu_lib_ok" == "xno" ; then 
    dnl Maybe lib?
    SLU_LIBS="$mld2p4_cv_superlu -L$mld2p4_cv_superludir/lib";
    LIBS="$SLU_LIBS -lm $save_LIBS";
    AC_TRY_LINK_FUNC(superlu_malloc, 
		     [mld2p4_cv_have_superlu=yes;pac_slu_lib_ok=yes;],
		     [mld2p4_cv_have_superlu=no;pac_slu_lib_ok=no; SLU_LIBS=""; SLU_INCLUDES=""])
 fi
 AC_MSG_RESULT($pac_slu_lib_ok)
fi
if test "x$pac_slu_header_ok" == "xyes" ; then 
   AC_MSG_CHECKING([for superlu version 5])
   AC_LANG_PUSH([C])
   AC_COMPILE_IFELSE(
       [AC_LANG_SOURCE([[#include "slu_ddefs.h"
			 int testdslu()
			 { SuperMatrix AC, *L, *U;
			   int *perm_r, *perm_c,  *etree,  panel_size, permc_spec, relax, info;
			   superlu_options_t options;   SuperLUStat_t stat;
			   GlobalLU_t Glu;   
			   dgstrf(&options, &AC, relax, panel_size, etree,
				  NULL, 0, perm_c, perm_r, L, U, &Glu, &stat, &info);               
			   
			 }]])],
       [ AC_MSG_RESULT([yes]);      pac_slu_version="5";],
       [ AC_MSG_RESULT([no]);      pac_slu_version="3_4";])
   AC_LANG_POP([C])
fi   

LIBS="$save_LIBS";
CPPFLAGS="$save_CPPFLAGS";
AC_LANG_POP([C])
])dnl 

dnl @synopsis PAC_CHECK_SUPERLU_Dist
dnl
dnl Will try to find the SUPERLU_Dist library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_SUPERLUDIST,
[AC_ARG_WITH(superludist, AC_HELP_STRING([--with-superludist=LIBNAME], [Specify the libname for SUPERLUDIST library. Requires you also specify SuperLU. Default: "-lsuperlu_dist"]),
        [mld2p4_cv_superludist=$withval],
        [mld2p4_cv_superludist='-lsuperlu_dist'])
AC_ARG_WITH(superludistdir, AC_HELP_STRING([--with-superludistdir=DIR], [Specify the directory for SUPERLUDIST library and includes.]),
        [mld2p4_cv_superludistdir=$withval],
        [mld2p4_cv_superludistdir=''])

AC_ARG_WITH(superludistincdir, AC_HELP_STRING([--with-superludistincdir=DIR], [Specify the directory for SUPERLUDIST includes.]),
        [mld2p4_cv_superludistincdir=$withval],
        [mld2p4_cv_superludistincdir=''])

AC_ARG_WITH(superludistlibdir, AC_HELP_STRING([--with-superludistlibdir=DIR], [Specify the directory for SUPERLUDIST library.]),
        [mld2p4_cv_superludistlibdir=$withval],
        [mld2p4_cv_superludistlibdir=''])

AC_LANG_PUSH([C])
save_LIBS="$LIBS"
save_CPPFLAGS="$CPPFLAGS"
save_CC="$CC"
CC=${MPICC}
CPP="${CC} -E"
if test "x$mld2p4_cv_superludistincdir" != "x"; then 
 AC_MSG_NOTICE([sludist dir $mld2p4_cv_superludistincdir]) 
 SLUDIST_INCLUDES="-I$mld2p4_cv_superludistincdir"
elif test "x$mld2p4_cv_superludistdir" != "x"; then 
 AC_MSG_NOTICE([sludist dir $mld2p4_cv_superludistdir]) 
 SLUDIST_INCLUDES="-I$mld2p4_cv_superludistdir"
fi
if test "x$mld2p4_cv_superludistlibdir" != "x"; then 
   SLUDIST_LIBS="-L$mld2p4_cv_superludistlibdir"
elif test "x$mld2p4_cv_superludistdir" != "x"; then 
   SLUDIST_LIBS="-L$mld2p4_cv_superludir"
fi

LIBS="$SLUDIST_LIBS $save_LIBS"
CPPFLAGS="$SLUDIST_INCLUDES $save_CPPFLAGS"

AC_CHECK_HEADERS([superlu_ddefs.h],
 [pac_sludist_header_ok=yes],
 [pac_sludist_header_ok=no; SLUDIST_INCLUDES=""])
if test "x$pac_sludist_header_ok" == "xno" ; then 
dnl Maybe Include or include subdirs? 
  unset ac_cv_header_superlu_ddefs_h
  SLUDIST_INCLUDES="-I$mld2p4_cv_superludistdir/include -I$mld2p4_cv_superludistdir/Include"
  CPPFLAGS="$SLUDIST_INCLUDES $save_CPPFLAGS"

 AC_CHECK_HEADERS([superlu_ddefs.h],
		 [pac_sludist_header_ok=yes],
		 [pac_sludist_header_ok=no; SLUDIST_INCLUDES=""])
fi

if test "x$pac_sludist_header_ok" == "xyes" ; then 
      SLUDIST_LIBS="$mld2p4_cv_superludist $SLUDIST_LIBS"
      LIBS="$SLUDIST_LIBS -lm $save_LIBS";
      AC_MSG_CHECKING([for superlu_malloc_dist in $SLUDIST_LIBS])
      AC_TRY_LINK_FUNC(superlu_malloc_dist, 
       [mld2p4_cv_have_superludist=yes;pac_sludist_lib_ok=yes;],
       [mld2p4_cv_have_superludist=no;pac_sludist_lib_ok=no; 
          SLUDIST_LIBS=""; ])
  if test "x$pac_sludist_lib_ok" == "xno" ; then 
     dnl Maybe lib?
     SLUDIST_LIBS="$mld2p4_cv_superludist  -L$mld2p4_cv_superludistdir/lib";
     LIBS="$SLUDIST_LIBS -lm $save_LIBS";
     AC_TRY_LINK_FUNC(superlu_malloc_dist, 
		     [mld2p4_cv_have_superludist=yes;pac_sludist_lib_ok=yes;],
		     [mld2p4_cv_have_superludist=no;pac_sludist_lib_ok=no; 
		      SLUDIST_LIBS="";SLUDIST_INCLUDES=""])
 fi
 AC_MSG_RESULT($pac_sludist_lib_ok)
 AC_MSG_CHECKING([for superlu_dist version 4])
 AC_LANG_PUSH([C])
 ac_cc=${MPICC-$CC}
 AC_COMPILE_IFELSE(
       [AC_LANG_SOURCE([[   #include "superlu_ddefs.h"
			    int testdslud()
			    {  LUstruct_t *LUstruct;
			       int n; 
			       LUstructInit(n, LUstruct);     
			    }]])],
       [ AC_MSG_RESULT([yes]);     pac_sludist_version="4";],
       [ AC_MSG_RESULT([no]);      pac_sludist_version="2_3";])
   AC_LANG_POP([C])
fi
 LIBS="$save_LIBS";
 CPPFLAGS="$save_CPPFLAGS";
 CC="$save_CC";
AC_LANG_POP([C])
])dnl 

dnl @synopsis PAC_CHECK_MUMPS
dnl
dnl Will try to find the MUMPS library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_MUMPS,
[AC_ARG_WITH(mumps, AC_HELP_STRING([--with-mumps=LIBNAME], [Specify the libname for MUMPS. Default: "-lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lpord"]),
        [mld2p4_cv_mumps=$withval],
        [mld2p4_cv_mumps='-lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lpord'])
 AC_ARG_WITH(mumpsdir, AC_HELP_STRING([--with-mumpsdir=DIR], [Specify the directory for MUMPS library and includes. Note: you will need to add auxiliary libraries with --extra-libs; this depends on how MUMPS was configured and installed, at a minimum you will need SCALAPACK and BLAS]),
        [mld2p4_cv_mumpsdir=$withval],
        [mld2p4_cv_mumpsdir=''])

AC_ARG_WITH(mumpsincdir, AC_HELP_STRING([--with-mumpsincdir=DIR], [Specify the directory for MUMPS includes.]),
        [mld2p4_cv_mumpsincdir=$withval],
        [mld2p4_cv_mumpsincdir=''])

AC_ARG_WITH(mumpslibdir, AC_HELP_STRING([--with-mumpslibdir=DIR], [Specify the directory for MUMPS library.]),
        [mld2p4_cv_mumpslibdir=$withval],
        [mld2p4_cv_mumpslibdir=''])

AC_LANG_PUSH([Fortran])
save_LIBS="$LIBS"
save_FC="$FC"
FC=${MPIFC}
if test "x$mld2p4_cv_mumpsincdir" != "x"; then 
 AC_MSG_NOTICE([mumps dir $mld2p4_cv_mumpsincdir]) 
 MUMPS_INCLUDES="$FMFLAG$mld2p4_cv_mumpsincdir"
elif test "x$mld2p4_cv_mumpsdir" != "x"; then 
 AC_MSG_NOTICE([mumps dir $mld2p4_cv_mumpsdir]) 
 MUMPS_INCLUDES="$FMFLAG$mld2p4_cv_mumpsdir"
fi
if test "x$mld2p4_cv_mumpslibdir" != "x"; then 
   MUMPS_LIBS="-L$mld2p4_cv_mumpslibdir"
elif test "x$mld2p4_cv_mumpsdir" != "x"; then 
   MUMPS_LIBS="-L$mld2p4_cv_mumpsdir"
fi

LIBS="$MUMPS_LIBS $save_LIBS $EXTRA_LIBS"
CPPFLAGS="$MUMPS_INCLUDES $save_CPPFLAGS"

ac_objext='o'
ac_ext='f90'
ac_fc="${MPIFC-$FC}";
save_FCFLAGS="$FCFLAGS";
FCFLAGS=" $MUMPS_INCLUDES $save_FCFLAGS"
AC_COMPILE_IFELSE([
		    program test
		    use dmumps_struc_def
		    end program test],
		  [pac_mumps_header_ok=yes; mld2p4_cv_mumpsincdir="$MUMPS_INCLUDES";],
		   [pac_mumps_header_ok=no; MUMPS_INCLUDES=""])
if test "x$pac_mumps_header_ok" == "xno" ; then 
   dnl Maybe Include or include subdirs? 
   MUMPS_INCLUDES="$FMFLAG$mld2p4_cv_mumpsdir/include"
   FCFLAGS="$MUMPS_INCLUDES $save_CPPFLAGS"
   
   AC_COMPILE_IFELSE([
		      program test
		      use dmumps_struc_def
		      end program test],
		     [pac_mumps_header_ok=yes mld2p4_cv_mumpsincdir="$MUMPS_INCLUDES";],
		     [pac_mumps_header_ok=no; MUMPS_INCLUDES=""])
   fi
if test "x$pac_mumps_header_ok" == "xno" ; then 
   dnl Maybe Include or include subdirs? 
   MUMPS_INCLUDES="$FMFLAG$mld2p4_cv_mumpsdir/Include"
   FCFLAGS="$MUMPS_INCLUDES $save_CPPFLAGS"
   
   AC_COMPILE_IFELSE([
		      program test
		      use dmumps_struc_def
		      end program test],
		     [pac_mumps_header_ok=yes mld2p4_cv_mumpsincdir="$MUMPS_INCLUDES";],
		     [pac_mumps_header_ok=no; MUMPS_INCLUDES=""])
   fi
   

if test "x$pac_mumps_header_ok" == "xyes" ; then 
      MUMPS_LIBS="$mld2p4_cv_mumps $MUMPS_LIBS"
      LIBS="$MUMPS_LIBS  $save_LIBS  $EXTRA_LIBS";
      AC_MSG_CHECKING([for dmumps in $MUMPS_LIBS])
      AC_TRY_LINK_FUNC(dmumps, 
       [mld2p4_cv_have_mumps=yes;pac_mumps_lib_ok=yes;],
       [mld2p4_cv_have_mumps=no;pac_mumps_lib_ok=no; 
          MUMPS_LIBS=""; ])
  if test "x$pac_mumps_lib_ok" == "xno" ; then 
     dnl Maybe lib?
     MUMPS_LIBS="$mld2p4_cv_mumps  -L$mld2p4_cv_mumpsdir/lib";
     LIBS="$MUMPS_LIBS  $save_LIBS  $EXTRA_LIBS";
     AC_TRY_LINK_FUNC(dmumps, 
		     [mld2p4_cv_have_mumps=yes;pac_mumps_lib_ok=yes;],
		     [mld2p4_cv_have_mumps=no;pac_mumps_lib_ok=no; 
		      MUMPS_LIBS="";MUMPS_INCLUDES=""])
 fi
 AC_MSG_RESULT($pac_mumps_lib_ok)
fi
 LIBS="$save_LIBS";
 CPPFLAGS="$save_CPPFLAGS";
 FC="$save_FC";
AC_LANG_POP([Fortran])
])dnl 



dnl @synopsis PAC_ARG_SERIAL_MPI
dnl
dnl Test for --with-serial-mpi={yes|no}
dnl 
dnl 
dnl
dnl Example use:
dnl
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_SERIAL_MPI],
[
AC_MSG_CHECKING([whether we want serial (fake) mpi])
AC_ARG_ENABLE(serial,
AC_HELP_STRING([--enable-serial], 
[Specify whether to enable a fake mpi library to run in serial mode. ]),
[
pac_cv_serial_mpi="yes";
]
dnl ,
dnl [pac_cv_serial_mpi="no";]
)
if test x"$pac_cv_serial_mpi" == x"yes" ; then
   AC_MSG_RESULT([yes.])
else
 pac_cv_serial_mpi="no";
 AC_MSG_RESULT([no.])
fi
]
)
