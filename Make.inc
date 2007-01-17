.mod=.mod
.fh=.fh
.SUFFIXES: .f90 $(.mod) .F90


#######################  Section 1 #######################
#  Define your compilers and compiler flags here         #
##########################################################
F90=/usr/local/gcc42/bin/gfortran
FC=/usr/local/gcc42/bin/gfortran
F77=$(FC)
CC=/usr/local/gcc42/bin/gcc
F90COPT=-O3 -ggdb
FCOPT=-O3 -ggdb
CCOPT=-O3 -ggdb

#######################  Section 2 #######################
#  Define your linker and linker flags here              #
##########################################################
F90LINK=/usr/local/mpich-gcc42/bin/mpif90 
FLINK=/usr/local/mpich-gcc42/bin/mpif77 
MPF90=/usr/local/mpich-gcc42/bin/mpif90 
MPF77=/usr/local/mpich-gcc42/bin/mpif77 	
MPCC=/usr/local/mpich-gcc42/bin/mpicc 

#######################  Section 3 #######################
#  Specify paths to libraries                            #
##########################################################
BLAS=-lblas-gcc42 -L$(HOME)/LIB
BLACS=-lmpiblacs-gcc42 -L$(HOME)/LIB


#######################  Section 4 #######################
#  Other useful tools&defines                            #
##########################################################
SLUDIR=/usr/local/SuperLU_3.0
SLU=-lslu_lx_gcc42 -L$(SLUDIR)
SLUDEF=-DHave_SLU_ -I$(SLUDIR)

UMFDIR=$(HOME)/LIB/Umfpack_gcc41
UMF=-lumfpack -lamd -L$(UMFDIR)
UMFDEF=-DHave_UMF_ -I$(UMFDIR)
#
#  We are using the public domain tool METIS from U. Minnesota. To get it
#  check URL http://www.cs.umn.edu:~karypis
#  
METIS_LIB = -L$(HOME)/NUMERICAL/metis-4.0 -lmetis 
LDLIBS=$(BLACS) $(SLU) $(UMF) $(BLAS) $(METIS_LIB)

#  Add -DLargeFptr for 64-bit addresses
CDEFINES=-DAdd_ $(SLUDEF) $(UMFDEF)
FDEFINES=-DNETLIB_BLACS -DHAVE_MOVE_ALLOC

AR=ar -cur
RANLIB=ranlib



#######################  Section 5 #######################
#  Do not edit this                                      #
##########################################################
LIBDIR=lib
BASELIBNAME=libpsb_base.a
PRECLIBNAME=libpsb_prec.a
METHDLIBNAME=libpsb_krylov.a
UTILLIBNAME=libpsb_util.a

# Under Linux/gmake there is a rule interpreting .mod as Modula source! 
$(.mod).o:
	
.f.o:
	$(FC) $(FCOPT) $(INCDIRS)  -c $<
.c.o:
	$(CC) $(CCOPT) $(INCDIRS) $(CDEFINES) -c $<
.f$(.mod):
	$(F90) $(FCOPT) $(INCDIRS) -c $<
.f90$(.mod):
	$(F90) $(F90COPT) $(INCDIRS) -c $<
.f90.o:
	$(F90) $(F90COPT) $(INCDIRS) -c $<
.F90.o:
	$(F90) $(F90COPT) $(INCDIRS) $(FDEFINES) -c $<



