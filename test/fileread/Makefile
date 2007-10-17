#PSBLAS_DIR=../psblas2-dev
#include $(PSBLAS_DIR)/Make.inc
#
# Libraries used
#
#LIBDIR=$(PSBLAS_DIR)/lib/
MLDDIR=../..
include $(MLDDIR)/Make.inc
PSBDIR=$(PSBLASDIR)/lib/
MLDLIBDIR=$(MLDDIR)/lib
MLD_LIB=-L$(MLDLIBDIR) -lmld_krylov -lmld_prec
PSBLAS_LIB= -L$(PSBDIR)  -lpsb_util  -lpsb_base

INCDIRS=-I$(MLDLIBDIR) -I$(PSBDIR)

DFOBJS=precdata.o getp.o df_bench.o enablecore.o
ZFOBJS=precdata.o getp.o zf_bench.o enablecore.o

EXEDIR=./runs

all: df_bench zf_bench

df_bench.o zf_bench.o: getp.o
getp.o: precdata.o 

df_bench: $(DFOBJS)
	$(F90LINK) $(LINKOPT) $(DFOBJS) -o df_bench \
	$(MLD_LIB)  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv df_bench $(EXEDIR) 

zf_bench: $(ZFOBJS)
	$(F90LINK) $(LINKOPT) $(ZFOBJS) -o zf_bench \
	$(MLD_LIB) $(PSBLAS_LIB)  $(LDLIBS)
	/bin/mv zf_bench $(EXEDIR) 


.f90.o:
	$(MPF90) $(F90COPT) $(INCDIRS) -c $<
#mmio.o: mmio.f90
#	$(MPF90) -O0        $(INCDIRS) -c $<


#$(DFLOBJS) $(ZFOBJS) $(DFOBJS):$(MODS)

clean: 
	/bin/rm -f $(DFOBJS) $(ZFOBJS)  \
        *$(.mod)  $(EXEDIR)/df_bench $(EXEDIR)/zf_bench

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

