MLDDIR=../..
include $(MLDDIR)/Make.inc
include $(MLDDIR)/Make.inc
PSBLIBDIR=$(PSBLASDIR)/lib/
PSBINCDIR=$(PSBLASDIR)/include
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIB=-L$(MLDLIBDIR) -lpsb_krylov -lmld_prec -lpsb_prec
PSBLAS_LIB= -L$(PSBLIBDIR)  -lpsb_util  -lpsb_base
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDLIBDIR)  $(FMFLAG)$(PSBINCDIR) $(FIFLAG).

PDOBJS=ppde.o  data_input.o mld_d_tlu_solver.o
PSOBJS=spde.o  data_input.o 
EXEDIR=./runs

all: ppde spde

ppde: $(PDOBJS)
	$(F90LINK) $(PDOBJS) -o ppde $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde $(EXEDIR)


spde: $(PSOBJS)
	$(F90LINK) $(PSOBJS) -o spde $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv spde $(EXEDIR)

ppde.o spde.o: data_input.o 
ppde.o: mld_d_tlu_solver.o


clean: 
	/bin/rm -f $(PDOBJS) $(EXEDIR)/ppde $(PSOBJS) $(EXEDIR)/spde 

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



