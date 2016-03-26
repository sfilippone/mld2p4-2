MLDDIR=../..
MLDINCDIR=$(MLDDIR)/include
include $(MLDINCDIR)/Make.inc.mld2p4
PSBLIBDIR=$(PSBLASDIR)/lib/
PSBINCDIR=$(PSBLASDIR)/include
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIB=-L$(MLDLIBDIR)  -lmld_prec
PSBLAS_LIB= -L$(PSBLIBDIR) -lpsb_krylov -lpsb_prec -lpsb_util  -lpsb_base
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDINCDIR)  $(FMFLAG)$(PSBINCDIR) $(FIFLAG).

DMOBJS=mld_dexample_ml.o data_input.o
D1OBJS=mld_dexample_1lev.o data_input.o
SMOBJS=mld_sexample_ml.o data_input.o
S1OBJS=mld_sexample_1lev.o data_input.o

EXEDIR=./runs

all: mld_sexample_ml mld_sexample_1lev mld_dexample_ml mld_dexample_1lev


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

mld_sexample_ml: $(SMOBJS)
	$(F90LINK) $(LINKOPT) $(SMOBJS) -o mld_sexample_ml \
	$(MLD_LIB) $(PSBLAS_LIB)  $(LDLIBS)
	/bin/mv mld_sexample_ml $(EXEDIR) 

mld_sexample_1lev: $(S1OBJS)
	$(F90LINK) $(LINKOPT) $(S1OBJS) -o mld_sexample_1lev \
	$(MLD_LIB) $(PSBLAS_LIB)  $(LDLIBS)
	/bin/mv mld_sexample_1lev $(EXEDIR) 

mld_sexample_ml.o: data_input.o
mld_sexample_1lev.o: data_input.o

.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) -c $<

clean: 
	/bin/rm -f $(DMOBJS) $(D1OBJS) $(SMOBJS) $(S1OBJS) \
        *$(.mod)  $(EXEDIR)/mld_dexample_ml $(EXEDIR)/mld_dexample_1lev\
	 $(EXEDIR)/mld_sexample_ml $(EXEDIR)/mld_sexample_1lev

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

