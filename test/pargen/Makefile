MLDDIR=../..
include $(MLDDIR)/Make.inc
PSBDIR=$(PSBLASDIR)/lib/
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIB=-L$(MLDLIBDIR) -lmld_krylov -lmld_prec
PSBLAS_LIB= -L$(PSBDIR)  -lpsb_util  -lpsb_base
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDLIBDIR)  $(FMFLAG)$(PSBDIR) $(FIFLAG).


EXEDIR=./runs


ppde: ppde.o  
	$(F90LINK) ppde.o  -o ppde $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde $(EXEDIR)


.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<


clean: 
	/bin/rm -f ppde.o   $(EXEDIR)/ppde 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



