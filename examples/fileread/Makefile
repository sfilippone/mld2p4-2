MLDDIR=../..
MLDINCDIR=$(MLDDIR)/include
include $(MLDINCDIR)/Make.inc.mld2p4
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIBS=-L$(MLDLIBDIR) -lpsb_krylov -lmld_prec -lpsb_prec
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDINCDIR) $(PSBLAS_INCLUDES) $(FIFLAG).

DMOBJS=mld_dexample_ml.o data_input.o
D1OBJS=mld_dexample_1lev.o data_input.o
ZMOBJS=mld_zexample_ml.o data_input.o
Z1OBJS=mld_zexample_1lev.o data_input.o
SMOBJS=mld_sexample_ml.o data_input.o
S1OBJS=mld_sexample_1lev.o data_input.o
CMOBJS=mld_cexample_ml.o data_input.o
C1OBJS=mld_cexample_1lev.o data_input.o

EXEDIR=./runs

all: mld_dexample_ml mld_dexample_1lev mld_zexample_ml mld_zexample_1lev\
	mld_sexample_ml mld_sexample_1lev mld_cexample_ml mld_cexample_1lev

mld_dexample_ml: $(DMOBJS)
	$(F90LINK) $(LINKOPT) $(DMOBJS) -o mld_dexample_ml \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_dexample_ml $(EXEDIR) 

mld_dexample_1lev: $(D1OBJS)
	$(F90LINK) $(LINKOPT) $(D1OBJS) -o mld_dexample_1lev \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_dexample_1lev $(EXEDIR) 

mld_dexample_ml.o: data_input.o
mld_dexample_1lev.o: data_input.o

mld_zexample_ml: $(ZMOBJS)
	$(F90LINK) $(LINKOPT) $(ZMOBJS) -o mld_zexample_ml \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_zexample_ml $(EXEDIR) 

mld_zexample_1lev: $(Z1OBJS)
	$(F90LINK) $(LINKOPT) $(Z1OBJS) -o mld_zexample_1lev \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_zexample_1lev $(EXEDIR) 

mld_zexample_ml.o: data_input.o
mld_zexample_1lev.o: data_input.o



mld_sexample_ml: $(SMOBJS)
	$(F90LINK) $(LINKOPT) $(SMOBJS) -o mld_sexample_ml \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_sexample_ml $(EXEDIR) 

mld_sexample_1lev: $(S1OBJS)
	$(F90LINK) $(LINKOPT) $(S1OBJS) -o mld_sexample_1lev \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_sexample_1lev $(EXEDIR) 

mld_sexample_ml.o: data_input.o
mld_sexample_1lev.o: data_input.o

mld_cexample_ml: $(CMOBJS)
	$(F90LINK) $(LINKOPT) $(CMOBJS) -o mld_cexample_ml \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_cexample_ml $(EXEDIR) 

mld_cexample_1lev: $(C1OBJS)
	$(F90LINK) $(LINKOPT) $(C1OBJS) -o mld_cexample_1lev \
	$(MLD_LIBS) $(PSBLAS_LIBS)  $(LDLIBS)
	/bin/mv mld_cexample_1lev $(EXEDIR) 

mld_cexample_ml.o: data_input.o
mld_cexample_1lev.o: data_input.o


.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) -c $<

clean: 
	/bin/rm -f *$(.mod) \
	$(DMOBJS) $(D1OBJS) $(ZMOBJS) $(Z1OBJS) \
        $(EXEDIR)/mld_dexample_ml $(EXEDIR)/mld_dexample_1lev\
	$(EXEDIR)/mld_zexample_ml $(EXEDIR)/mld_zexample_1lev\
	$(SMOBJS) $(S1OBJS) $(CMOBJS) $(C1OBJS) \
        $(EXEDIR)/mld_sexample_ml $(EXEDIR)/mld_sexample_1lev\
	$(EXEDIR)/mld_cexample_ml $(EXEDIR)/mld_cexample_1lev

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

