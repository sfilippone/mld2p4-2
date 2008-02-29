include ../Make.inc
LIBDIR=../lib
PSBLIBDIR=$(PSBLASDIR)/lib
HERE=.

FINCLUDES=$(FMFLAG). $(FMFLAG)$(LIBDIR) $(FMFLAG)$(PSBLIBDIR)


MODOBJS=mld_prec_type.o  mld_prec_mod.o mld_inner_mod.o mld_basep_bld_mod.o
MPFOBJS=mld_daggrmat_raw_asb.o mld_daggrmat_smth_asb.o \
	mld_zaggrmat_raw_asb.o mld_zaggrmat_smth_asb.o 
MPCOBJS=mld_dslud_interface.o mld_zslud_interface.o
INNEROBJS=mld_das_bld.o mld_dslu_bld.o mld_dumf_bld.o mld_dilu0_fact.o\
	mld_dmlprec_bld.o mld_dsp_renum.o mld_dfact_bld.o mld_dilu_bld.o \
	mld_dbaseprec_bld.o mld_ddiag_bld.o mld_daggrmap_bld.o \
	mld_dmlprec_aply.o mld_dslud_bld.o\
	mld_dbaseprec_aply.o mld_dsub_aply.o mld_dsub_solve.o \
	mld_das_aply.o mld_daggrmat_asb.o \
	mld_zas_bld.o mld_zslu_bld.o mld_zumf_bld.o mld_zilu0_fact.o\
	mld_zmlprec_bld.o mld_zsp_renum.o mld_zfact_bld.o mld_zilu_bld.o \
	mld_zbaseprec_bld.o mld_zdiag_bld.o mld_zaggrmap_bld.o \
	mld_zmlprec_aply.o  mld_zslud_bld.o\
	mld_zbaseprec_aply.o mld_zsub_aply.o mld_zsub_solve.o \
	mld_zas_aply.o mld_zaggrmat_asb.o\
	mld_diluk_fact.o mld_ziluk_fact.o mld_dilut_fact.o mld_zilut_fact.o \
	$(MPFOBJS) 
OUTEROBJS=mld_dprecbld.o mld_dprecfree.o mld_dprecset.o mld_dprecinit.o\
	mld_dprec_aply.o \
	mld_zprecbld.o mld_zprecfree.o mld_zprecset.o mld_zprecinit.o \
	mld_zprec_aply.o 
F90OBJS=$(OUTEROBJS) $(INNEROBJS)

COBJS=mld_dslu_interface.o mld_dumf_interface.o mld_zslu_interface.o mld_zumf_interface.o
OBJS=$(F90OBJS) $(COBJS) $(MPFOBJS) $(MPCOBJS) $(MODOBJS)

LIBMOD=mld_prec_mod$(.mod) 
LOCAL_MODS=$(LIBMOD) mld_prec_type$(.mod) mld_inner_mod$(.mod) 
LIBNAME=libmld_prec.a

lib: mpobjs $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(LIBMOD) $(LOCAL_MODS) $(LIBDIR)

$(F90OBJS) $(MPFOBJS): $(MODOBJS:.o=$(.mod))
mld_prec_mod.o mld_innner_mod.o: mld_prec_type.o
mld_inner_mod.o: mld_basep_bld_mod.o

$(MODOBJS): $(PSBLIBDIR)/psb_base_mod$(.mod)
 
mpobjs: $(MODOBJS)
	(make $(MPFOBJS) F90="$(MPF90)" F90COPT="$(F90COPT)")
	(make $(MPCOBJS) CC="$(MPCC)" CCOPT="$(CCOPT)")

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)
