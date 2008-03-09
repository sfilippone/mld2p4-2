##########################################################
#                                                        #
#  Note: directories external to the MLD2P4 subtree      #
#  must be specified here with absolute pathnames        #
#                                                        #
##########################################################
PSBLASDIR=/home/sfilippo/NUMERICAL/PSBLAS2/psblas2-dev/
include $(PSBLASDIR)/Make.inc


INSTALL=/usr/bin/install -c
INSTALL_DATA=${INSTALL} -m 644
INSTALL_DIR=/usr/local/PSBLAS/MLD2P4/gnu43
INSTALL_LIBDIR=/usr/local/PSBLAS/MLD2P4/gnu43/lib
INSTALL_INCLUDEDIR=/usr/local/PSBLAS/MLD2P4/gnu43/include

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

SLULIBS=
SLUFLAGS=

SLUDISTLIBS=
SLUDISTFLAGS=


UMFLIBS=-lumfpack -lamd -L/home/sfilippo/LIB/Umfpack_gcc41/
UMFFLAGS=-DHave_UMF_ -I/home/sfilippo/LIB/Umfpack_gcc41/


LDLIBS=$(SLULIBS) $(SLUDISTLIBS) $(UMFLIBS) $(PSBLDLIBS)

#
CDEFINES=$(SLUFLAGS) $(UMFFLAGS) $(SLUDISTFLAGS) $(PSBCDEFINES)
FDEFINES=$(PSBFDEFINES)
