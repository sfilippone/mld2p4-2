MLDDIR=../..
include $(MLDDIR)/Make.inc
PSBDIR=$(PSBLASDIR)/lib/
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIB=-L$(MLDLIBDIR) -lmld_krylov -lmld_prec
PSBLAS_LIB= -L$(PSBDIR)  -lpsb_util  -lpsb_base
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDLIBDIR)  $(FMFLAG)$(PSBDIR) $(FIFLAG).

DFOBJS=df_bench.o 
DFSOBJS=df_sample.o data_input.o
SFSOBJS=sf_sample.o data_input.o
CFSOBJS=cf_sample.o data_input.o
ZFSOBJS=zf_sample.o data_input.o
ZFOBJS=zf_bench.o 

EXEDIR=./runs

all: df_bench zf_bench sf_sample df_sample cf_sample zf_sample

df_bench: $(DFOBJS)
	$(F90LINK) $(LINKOPT) $(DFOBJS) -o df_bench \
	$(MLD_LIB)  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv df_bench $(EXEDIR) 

df_sample: $(DFSOBJS)
	$(F90LINK) $(LINKOPT) $(DFSOBJS) -o df_sample \
	$(MLD_LIB)  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv df_sample $(EXEDIR) 

sf_sample: $(SFSOBJS)
	$(F90LINK) $(LINKOPT) $(SFSOBJS) -o sf_sample \
	$(MLD_LIB)  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv sf_sample $(EXEDIR) 

cf_sample: $(CFSOBJS)
	$(F90LINK) $(LINKOPT) $(CFSOBJS) -o cf_sample \
	$(MLD_LIB)  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv cf_sample $(EXEDIR) 

zf_sample: $(ZFSOBJS)
	$(F90LINK) $(LINKOPT) $(ZFSOBJS) -o zf_sample \
	$(MLD_LIB)  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv zf_sample $(EXEDIR) 

zf_bench: $(ZFOBJS)
	$(F90LINK) $(LINKOPT) $(ZFOBJS) -o zf_bench \
	$(MLD_LIB) $(PSBLAS_LIB)  $(LDLIBS)
	/bin/mv zf_bench $(EXEDIR) 

sf_sample.o: data_input.o
df_sample.o: data_input.o
cf_sample.o: data_input.o
zf_sample.o: data_input.o

.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) -c $<

clean: 
	/bin/rm -f $(DFOBJS) $(ZFOBJS) $(SFSOBJS)  $(DFSOBJS) \
        *$(.mod)  $(EXEDIR)/df_bench $(EXEDIR)/zf_bench \
	$(EXEDIR)/sf_sample $(EXEDIR)/cf_sample  \
	$(EXEDIR)/df_sample $(EXEDIR)/zf_sample

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

