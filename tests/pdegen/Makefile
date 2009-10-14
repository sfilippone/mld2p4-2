MLDDIR=../..
include $(MLDDIR)/Make.inc
PSBDIR=$(PSBLASDIR)/lib/
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIB=-L$(MLDLIBDIR)  -lmld_prec
PSBLAS_LIB= -L$(PSBDIR) -lpsb_krylov -lpsb_prec  -lpsb_util  -lpsb_base
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDLIBDIR)  $(FMFLAG)$(PSBDIR) $(FIFLAG).


EXEDIR=./runs

all: ppde spde

ppde: ppde.o  data_input.o
	$(F90LINK) ppde.o  data_input.o -o ppde $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde $(EXEDIR)

spde: spde.o   data_input.o
	$(F90LINK) spde.o  data_input.o -o spde $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv spde $(EXEDIR)

ppde.o spde.o: data_input.o 

.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<

clean: 
	/bin/rm -f data_input.o ppde.o $(EXEDIR)/ppde  spde.o $(EXEDIR)/spde 

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



