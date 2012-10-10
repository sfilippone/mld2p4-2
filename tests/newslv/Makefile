MLDDIR=../..
include $(MLDDIR)/Make.inc
PSBLIBDIR=$(PSBLASDIR)/lib/
PSBINCDIR=$(PSBLASDIR)/include
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIB=-L$(MLDLIBDIR) -lpsb_krylov -lmld_prec -lpsb_prec
PSBLAS_LIB= -L$(PSBLIBDIR)  -lpsb_util  -lpsb_base
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDLIBDIR)  $(FMFLAG)$(PSBINCDIR) $(FIFLAG).

PD3DOBJS=ppde3d.o  data_input.o mld_d_tlu_solver.o mld_d_tlu_solver_impl.o
PSOBJS=spde.o  data_input.o 
EXEDIR=./runs

all: ppde3d 


ppde3d: $(PD3DOBJS)
	$(F90LINK) $(PD3DOBJS) -o ppde3d $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde3d $(EXEDIR)


ppde3d.o mld_d_tlu_solver_impl.o: data_input.o mld_d_tlu_solver.o


clean: 
	/bin/rm -f $(PD3DOBJS) $(EXEDIR)/ppde3d 

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



