##########################################################
#                                                        #
#  Note: directories external to the MLD2P4 subtree      #
#  must be specified here with absolute pathnames        #
#                                                        #
##########################################################
PSBLASDIR=$(HOME)/NUMERICAL/PSBLAS2/psblas2-dev
include $(PSBLASDIR)/Make.inc

##########################################################
#                                                        #
#  Additional defines and libraries for multilevel       #
#  Note that these libraries should be compatible        #
#  (compiled with) the compilers specified in the        #
#  PSBLAS main Make.inc                                  #
#                                                        #
#                                                        #
##########################################################
SLUINCDIR=/usr/local/SuperLU_3.0
SLULIBDIR=/usr/local/SuperLU_3.0
SLU=-lslu_lx_gcc42 -L$(SLULIBDIR)
#SLUDistLIBDIR=/usr/local/SLUDist_2.0
#SLUDistINCDIR=/usr/local/SLUDist_2.0
#SLUDist=-lslud-gcc42 -L$(SLUDistLIBDIR) -L $(HOME)/LIB
#SLUDISTDEF=-DHave_SLUDist_ -I$(SLUDistINCDIR)
SLUDEF=-DHave_SLU_ -I$(SLUINCDIR) 


UMFINCDIR=$(HOME)/LIB/Umfpack_gcc41
UMFLIBDIR=$(HOME)/LIB/Umfpack_gcc41
UMF=-lumfpack -lamd -L$(UMFLIBDIR)
UMFDEF=-DHave_UMF_ -I$(UMFINCDIR)

LDLIBS=$(BLACS) $(SLU) $(SLUDist) $(UMF) $(BLAS) $(METIS_LIB)

#  Add -DLargeFptr for 64-bit addresses
CDEFINES=-DAdd_ $(SLUDEF) $(UMFDEF) $(SLUDISTDEF)
