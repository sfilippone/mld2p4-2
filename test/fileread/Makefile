MLDDIR=../..
include $(MLDDIR)/Make.inc
PSBDIR=$(PSBLASDIR)/lib/
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIB=-L$(MLDLIBDIR) -lmld_krylov -lmld_prec
PSBLAS_LIB= -L$(PSBDIR)  -lpsb_util  -lpsb_base
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MLDLIBDIR)  $(FMFLAG)$(PSBDIR) $(FIFLAG).

DFOBJS=df_bench.o 
DFSOBJS=df_sample.o
ZFOBJS=zf_bench.o 

EXEDIR=./runs

all: df_bench zf_bench df_sample

df_bench: $(DFOBJS)
	$(F90LINK) $(LINKOPT) $(DFOBJS) -o df_bench \
	$(MLD_LIB)  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv df_bench $(EXEDIR) 

df_sample: $(DFSOBJS)
	$(F90LINK) $(LINKOPT) $(DFSOBJS) -o df_sample \
	$(MLD_LIB)  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv df_sample $(EXEDIR) 

zf_bench: $(ZFOBJS)
	$(F90LINK) $(LINKOPT) $(ZFOBJS) -o zf_bench \
	$(MLD_LIB) $(PSBLAS_LIB)  $(LDLIBS)
	/bin/mv zf_bench $(EXEDIR) 


.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) -c $<

clean: 
	/bin/rm -f $(DFOBJS) $(ZFOBJS)  \
        *$(.mod)  $(EXEDIR)/df_bench $(EXEDIR)/zf_bench

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

