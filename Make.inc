##########################################################
#                                                        #
#  Note: directories external to the MLD2P4 subtree      #
#  must be specified here with absolute pathnames        #
#                                                        #
##########################################################
PSBLASDIR=/home/sfilippo/NUMERICAL/PSBLAS2/psblas-2.2-maint
include $(PSBLASDIR)/Make.inc

##########################################################
#                                                        #
#  Additional defines and libraries for multilevel       #
#  Note that these libraries should be compatible        #
#  (compiled with) the compilers specified in the        #
#  PSBLAS main Make.inc                                  #
#                                                        #
#  Examples:                                             #
#  UMFLIBS=-lumfpack -lamd -L/path/to/UMFPACK            #
#  UMFFLAGS=-DHave_UMF_ -I/path/to/UMFPACK               #
#                                                        #
#  SLULIBS=-lslu   -L/path/to/SuperLU                    #
#  SLUFLAGS=-DHave_SLU_ -I/path/to/SuperLU               #
#                                                        #
#  SLUDISTLIBS=-lslud -L/path/to/SuperLUDist             #
#  SLUDISTFLAGS=-DHave_SLUDist_ -I/path/to/SuperLUDist   #
#                                                        #
##########################################################

SLULIBS=-lslu -L/usr/local/SuperLU_3.0/gnu43
SLUFLAGS=-DHave_SLU_ -I/usr/local/SuperLU_3.0/gnu43

SLUDISTLIBS=-lslud -L/usr/local/SLUDist_2.0/gnu43
SLUDISTFLAGS=-DHave_SLUDist_ -I/usr/local/SLUDist_2.0/gnu43


UMFLIBS=-lumfpack -lamd -L/home/sfilippo/LIB/Umfpack_gcc41/
UMFFLAGS=-DHave_UMF_ -I/home/sfilippo/LIB/Umfpack_gcc41/


LDLIBS=$(SLULIBS) $(SLUDISTLIBS) $(UMFLIBS) $(PSBLDLIBS)

#
CDEFINES=$(SLUFLAGS) $(UMFFLAGS) $(SLUDISTFLAGS) $(PSBCDEFINES)
FDEFINES=$(PSBFDEFINES)
