include ../Make.inc
LIBDIR=../lib
PSBLIBDIR=$(PSBLASDIR)/lib
HERE=.

FINCLUDES=$(FMFLAG). $(FMFLAG)$(LIBDIR) $(FMFLAG)$(PSBLIBDIR)


MODOBJS=mld_prec_type.o  mld_prec_mod.o mld_inner_mod.o mld_basep_bld_mod.o
MPFOBJS=mld_saggrmat_raw_asb.o mld_saggrmat_smth_asb.o \
	mld_daggrmat_raw_asb.o mld_daggrmat_smth_asb.o \
	mld_caggrmat_raw_asb.o mld_caggrmat_smth_asb.o \
	mld_zaggrmat_raw_asb.o mld_zaggrmat_smth_asb.o 
MPCOBJS=mld_sslud_interface.o mld_dslud_interface.o mld_cslud_interface.o mld_zslud_interface.o
INNEROBJS=mld_sas_bld.o mld_sslu_bld.o mld_sumf_bld.o mld_silu0_fact.o\
	mld_smlprec_bld.o mld_ssp_renum.o mld_sfact_bld.o mld_silu_bld.o \
	mld_sbaseprec_bld.o mld_sdiag_bld.o mld_saggrmap_bld.o \
	mld_smlprec_aply.o mld_sslud_bld.o\
	mld_sbaseprec_aply.o mld_ssub_aply.o mld_ssub_solve.o \
	mld_sas_aply.o mld_saggrmat_asb.o \
	mld_das_bld.o mld_dslu_bld.o mld_dumf_bld.o mld_dilu0_fact.o\
	mld_dmlprec_bld.o mld_dsp_renum.o mld_dfact_bld.o mld_dilu_bld.o \
	mld_dbaseprec_bld.o mld_ddiag_bld.o mld_daggrmap_bld.o \
	mld_dmlprec_aply.o mld_dslud_bld.o\
	mld_dbaseprec_aply.o mld_dsub_aply.o mld_dsub_solve.o \
	mld_das_aply.o mld_daggrmat_asb.o \
	mld_cas_bld.o mld_cslu_bld.o mld_cumf_bld.o mld_cilu0_fact.o\
	mld_cmlprec_bld.o mld_csp_renum.o mld_cfact_bld.o mld_cilu_bld.o \
	mld_cbaseprec_bld.o mld_cdiag_bld.o mld_caggrmap_bld.o \
	mld_cmlprec_aply.o  mld_cslud_bld.o\
	mld_cbaseprec_aply.o mld_csub_aply.o mld_csub_solve.o \
	mld_cas_aply.o mld_caggrmat_asb.o\
	mld_zas_bld.o mld_zslu_bld.o mld_zumf_bld.o mld_zilu0_fact.o\
	mld_zmlprec_bld.o mld_zsp_renum.o mld_zfact_bld.o mld_zilu_bld.o \
	mld_zbaseprec_bld.o mld_zdiag_bld.o mld_zaggrmap_bld.o \
	mld_zmlprec_aply.o  mld_zslud_bld.o\
	mld_zbaseprec_aply.o mld_zsub_aply.o mld_zsub_solve.o \
	mld_zas_aply.o mld_zaggrmat_asb.o\
	mld_siluk_fact.o mld_ciluk_fact.o mld_silut_fact.o mld_cilut_fact.o \
	mld_diluk_fact.o mld_ziluk_fact.o mld_dilut_fact.o mld_zilut_fact.o \
	$(MPFOBJS) 
OUTEROBJS=mld_sprecbld.o mld_sprecfree.o mld_sprecset.o mld_sprecinit.o\
	mld_sprec_aply.o \
	mld_dprecbld.o mld_dprecfree.o mld_dprecset.o mld_dprecinit.o\
	mld_dprec_aply.o \
	mld_cprecbld.o mld_cprecfree.o mld_cprecset.o mld_cprecinit.o\
	mld_cprec_aply.o \
	mld_zprecbld.o mld_zprecfree.o mld_zprecset.o mld_zprecinit.o \
	mld_zprec_aply.o 
F90OBJS=$(OUTEROBJS) $(INNEROBJS)

COBJS= mld_sslu_interface.o mld_sumf_interface.o \
       mld_dslu_interface.o mld_dumf_interface.o \
       mld_cslu_interface.o mld_cumf_interface.o \
       mld_zslu_interface.o mld_zumf_interface.o
OBJS=$(F90OBJS) $(COBJS) $(MPFOBJS) $(MPCOBJS) $(MODOBJS)

LIBMOD=mld_prec_mod$(.mod) 
LOCAL_MODS=$(LIBMOD) mld_prec_type$(.mod) mld_inner_mod$(.mod) mld_basep_bld_mod$(.mod)
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
