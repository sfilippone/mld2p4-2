MLDDIR=../..
MLDINCDIR=$(MLDDIR)/include
include $(MLDINCDIR)/Make.inc.mld2p4
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIBS=-L$(MLDLIBDIR) -lpsb_krylov -lmld_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).

PD3DOBJS=ppde3d.o  data_input.o mld_d_tlu_solver.o mld_d_tlu_solver_impl.o
PSOBJS=spde.o  data_input.o 
EXEDIR=./runs

all: ppde3d 


ppde3d: $(PD3DOBJS)
	$(F90LINK) $(PD3DOBJS) -o ppde3d $(MLD_LIBS) $(PSBLAS_LIBS) $(LDLIBS)
	/bin/mv ppde3d $(EXEDIR)


ppde3d.o mld_d_tlu_solver_impl.o: data_input.o mld_d_tlu_solver.o


clean: 
	/bin/rm -f $(PD3DOBJS) $(EXEDIR)/ppde3d 

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



