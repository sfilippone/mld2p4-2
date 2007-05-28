PSBLASDIR=../psblas2
include $(PSBLASDIR)/Make.inc


LIBDIR=$(PSBLASDIR)/lib
HERE=.
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).


MODOBJS=psb_prec_type.o  psb_prec_mod.o
MPFOBJS=mld_daggrmat_asb.o mld_zaggrmat_asb.o 
MPCOBJS=mld_slud_impl.o mld_zslud_impl.o
F90OBJS=mld_dasmat_bld.o mld_dslu_bld.o mld_dumf_bld.o mld_dilu_fct.o\
	mld_dmlprec_bld.o mld_dsp_renum.o mld_dbjac_bld.o mld_dilu_bld.o \
	psb_dprecbld.o psb_dprecfree.o psb_dprecset.o psb_dprecinit.o\
	mld_dbaseprec_bld.o mld_ddiagsc_bld.o mld_daggrmap_bld.o \
	mld_dprec_aply.o mld_dmlprec_aply.o mld_dslud_bld.o\
	mld_dbaseprec_aply.o mld_dbjac_aply.o\
	mld_zasmat_bld.o mld_zslu_bld.o mld_zumf_bld.o mld_zilu_fct.o\
	mld_zmlprec_bld.o mld_zsp_renum.o mld_zbjac_bld.o mld_zilu_bld.o \
	psb_zprecbld.o psb_zprecfree.o psb_zprecset.o psb_zprecinit.o \
	mld_zbaseprec_bld.o mld_zdiagsc_bld.o mld_zaggrmap_bld.o \
	mld_zprec_aply.o mld_zmlprec_aply.o  mld_zslud_bld.o\
	mld_zbaseprec_aply.o mld_zbjac_aply.o\
	$(MPFOBJS) 
COBJS=mld_slu_impl.o mld_umf_impl.o mld_zslu_impl.o mld_zumf_impl.o
OBJS=$(F90OBJS) $(COBJS) $(MPFOBJS) $(MPCOBJS) $(MODOBJS)

LIBMOD=psb_prec_mod$(.mod)
LOCAL_MODS=$(LIBMOD) psb_prec_type$(.mod)
LIBNAME=$(PRECLIBNAME)

lib: mpobjs $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(LIBMOD) $(LOCAL_MODS) $(LIBDIR)

$(F90OBJS) $(MPFOBJS): $(MODOBJS)
psb_prec_mod.o: psb_prec_type.o
 
mpobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" F90COPT="$(F90COPT)")
	(make $(MPCOBJS) CC="$(MPCC)" CCOPT="$(CCOPT)")

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)
