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

all: mld_s_pde3d mld_d_pde3d mld_s_pde2d mld_d_pde2d 

mld_d_pde3d: mld_d_pde3d.o  data_input.o
	$(F90LINK) mld_d_pde3d.o  data_input.o -o mld_d_pde3d $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv mld_d_pde3d $(EXEDIR)

mld_s_pde3d: mld_s_pde3d.o   data_input.o
	$(F90LINK) mld_s_pde3d.o  data_input.o -o mld_s_pde3d $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv mld_s_pde3d $(EXEDIR)

mld_d_pde2d: mld_d_pde2d.o  data_input.o
	$(F90LINK) -g mld_d_pde2d.o  data_input.o -o mld_d_pde2d $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv mld_d_pde2d $(EXEDIR)

mld_s_pde2d: mld_s_pde2d.o   data_input.o
	$(F90LINK) mld_s_pde2d.o  data_input.o -o mld_s_pde2d $(MLD_LIB) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv mld_s_pde2d $(EXEDIR)

mld_d_pde3d.o mld_s_pde3d.o mld_d_pde2d.o mld_s_pde2d.o: data_input.o 

check: all
	cd runs && ./mld_d_pde2d <mld_d_pde.inp && ./mld_s_pde2d<mld_d_pde.inp


clean: 
	/bin/rm -f data_input.o mld_d_pde3d.o mld_s_pde3d.o mld_d_pde2d.o mld_s_pde2d.o *$(.mod)\
        $(EXEDIR)/mld_d_pde3d  $(EXEDIR)/mld_s_pde3d $(EXEDIR)/mld_d_pde2d  $(EXEDIR)/mld_s_pde2d 

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



