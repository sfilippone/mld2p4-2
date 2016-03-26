MLDDIR=../..
MLDINCDIR=$(MLDDIR)/include
include $(MLDINCDIR)/Make.inc.mld2p4
PSBLIBDIR=$(PSBLASDIR)/lib/
PSBINCDIR=$(PSBLASDIR)/include
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIB=-L$(MLDLIBDIR) -lpsb_krylov -lmld_prec -lpsb_prec
PSBLAS_LIB= -L$(PSBLIBDIR)  -lpsb_util  -lpsb_base
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDINCDIR)  $(FMFLAG)$(PSBINCDIR) $(FIFLAG).


EXEDIR=./runs

all: spde3d ppde3d spde2d ppde2d 

ppde3d: ppde3d.o  data_input.o
	$(F90LINK) ppde3d.o  data_input.o -o ppde3d $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde3d $(EXEDIR)

spde3d: spde3d.o   data_input.o
	$(F90LINK) spde3d.o  data_input.o -o spde3d $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv spde3d $(EXEDIR)


ppde2d: ppde2d.o  data_input.o
	$(F90LINK) -g ppde2d.o  data_input.o -o ppde2d $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde2d $(EXEDIR)

spde2d: spde2d.o   data_input.o
	$(F90LINK) spde2d.o  data_input.o -o spde2d $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv spde2d $(EXEDIR)

ppde3d.o spde3d.o ppde2d.o spde2d.o: data_input.o 

check: all
	cd runs && ./ppde2d <ppde.inp && ./spde2d<ppde.inp


clean: 
	/bin/rm -f data_input.o ppde3d.o spde3d.o ppde2d.o spde2d.o *$(.mod)\
        $(EXEDIR)/ppde3d  $(EXEDIR)/spde3d $(EXEDIR)/ppde2d  $(EXEDIR)/spde2d 

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



