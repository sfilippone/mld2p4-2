include ../../../Make.inc
LIBDIR=../../../lib
INCDIR=../../../include
PSBLIBDIR=$(PSBLASDIR)/lib
PSBINCDIR=$(PSBLASDIR)/include
HERE=../..

FINCLUDES=$(FMFLAG)../.. $(FMFLAG)$(INCDIR) $(FMFLAG)$(PSBINCDIR) $(FMFLAG)$(PSBLIBDIR)


OBJS=mld_d_base_aggregator_mat_asb.o \
mld_d_base_aggregator_tprol.o \
mld_d_hybrid_aggregator_tprol.o \
mld_d_symdec_aggregator_tprol.o  \
mld_d_map_to_tprol.o mld_d_dec_map_bld.o mld_d_hyb_map_bld.o

LIBNAME=libmld_prec.a

lib:  $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)

mpobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" F90COPT="$(F90COPT)")
	(make $(MPCOBJS) CC="$(MPCC)" CCOPT="$(CCOPT)")

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)
