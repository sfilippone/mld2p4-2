MLDDIR=../..
include $(MLDDIR)/Make.inc
PSBDIR=$(PSBLASDIR)/lib/
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIB=-L$(MLDLIBDIR) -lmld_krylov -lmld_prec
PSBLAS_LIB= -L$(PSBDIR)  -lpsb_util  -lpsb_base
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDLIBDIR)  $(FMFLAG)$(PSBDIR) $(FIFLAG).

DMOBJS=mld_dexample_ml.o data_input.o
D1OBJS=mld_dexample_1lev.o data_input.o

EXEDIR=./runs

all: mld_dexample_ml mld_dexample_1lev

mld_dexample_ml: $(DMOBJS)
	$(F90LINK) $(LINKOPT) $(DMOBJS) -o mld_dexample_ml \
	$(MLD_LIB) $(PSBLAS_LIB)  $(LDLIBS)
	/bin/mv mld_dexample_ml $(EXEDIR) 

mld_dexample_1lev: $(D1OBJS)
	$(F90LINK) $(LINKOPT) $(D1OBJS) -o mld_dexample_1lev \
	$(MLD_LIB) $(PSBLAS_LIB)  $(LDLIBS)
	/bin/mv mld_dexample_1lev $(EXEDIR) 

mld_dexample_ml.o: data_input.o
mld_dexample_1lev.o: data_input.o

.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) -c $<

clean: 
	/bin/rm -f $(DMOBJS) $(D1OBJS) \
        *$(.mod)  $(EXEDIR)/mld_dexample_ml $(EXEDIR)/mld_dexample_1lev

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

