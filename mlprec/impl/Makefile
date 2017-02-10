include ../../Make.inc
LIBDIR=../../lib
INCDIR=../../include
HERE=..

FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(INCDIR) $(PSBLAS_INCLUDES)



DMPFOBJS=mld_daggrmat_nosmth_asb.o mld_daggrmat_smth_asb.o mld_daggrmat_minnrg_asb.o mld_daggrmat_biz_asb.o

SMPFOBJS=mld_saggrmat_nosmth_asb.o mld_saggrmat_smth_asb.o mld_saggrmat_minnrg_asb.o mld_saggrmat_biz_asb.o

ZMPFOBJS=mld_zaggrmat_nosmth_asb.o mld_zaggrmat_smth_asb.o mld_zaggrmat_minnrg_asb.o mld_zaggrmat_biz_asb.o

CMPFOBJS=mld_caggrmat_nosmth_asb.o mld_caggrmat_smth_asb.o mld_caggrmat_minnrg_asb.o mld_caggrmat_biz_asb.o


MPFOBJS=$(SMPFOBJS) $(DMPFOBJS) $(CMPFOBJS) $(ZMPFOBJS) 

MPCOBJS=mld_sslud_interface.o mld_dslud_interface.o mld_cslud_interface.o mld_zslud_interface.o


DINNEROBJS= mld_dmlprec_bld.o \
	mld_d_smoothers_bld.o mld_d_hierarchy_bld.o \
	mld_dilu0_fact.o   mld_diluk_fact.o  mld_dilut_fact.o  mld_daggrmap_bld.o \
	mld_d_dec_map_bld.o mld_dmlprec_aply.o mld_daggrmat_asb.o   \
	$(DMPFOBJS) mld_d_extprol_bld.o mld_d_lev_aggrmap_bld.o  mld_d_lev_aggrmat_asb.o

SINNEROBJS=   mld_smlprec_bld.o \
	mld_s_smoothers_bld.o mld_s_hierarchy_bld.o \
	mld_silu0_fact.o   mld_siluk_fact.o  mld_silut_fact.o  mld_saggrmap_bld.o \
	mld_s_dec_map_bld.o mld_smlprec_aply.o mld_saggrmat_asb.o   \
	$(SMPFOBJS) mld_s_extprol_bld.o mld_s_lev_aggrmap_bld.o  mld_s_lev_aggrmat_asb.o

ZINNEROBJS=   mld_zmlprec_bld.o \
	mld_z_smoothers_bld.o mld_z_hierarchy_bld.o \
	mld_zilu0_fact.o   mld_ziluk_fact.o  mld_zilut_fact.o  mld_zaggrmap_bld.o \
	mld_z_dec_map_bld.o mld_zmlprec_aply.o mld_zaggrmat_asb.o   \
	$(ZMPFOBJS)  mld_z_extprol_bld.o mld_z_lev_aggrmap_bld.o  mld_z_lev_aggrmat_asb.o

CINNEROBJS=   mld_cmlprec_bld.o \
	mld_c_smoothers_bld.o mld_c_hierarchy_bld.o \
	mld_cilu0_fact.o   mld_ciluk_fact.o  mld_cilut_fact.o  mld_caggrmap_bld.o \
	mld_c_dec_map_bld.o mld_cmlprec_aply.o mld_caggrmat_asb.o   \
	$(CMPFOBJS)  mld_c_extprol_bld.o mld_c_lev_aggrmap_bld.o  mld_c_lev_aggrmat_asb.o

INNEROBJS= $(SINNEROBJS) $(DINNEROBJS) $(CINNEROBJS) $(ZINNEROBJS)


DOUTEROBJS=mld_dprecbld.o  mld_dprecset.o mld_dprecinit.o mld_dprecaply.o mld_dcprecset.o

SOUTEROBJS=mld_sprecbld.o  mld_sprecset.o mld_sprecinit.o mld_sprecaply.o mld_scprecset.o

ZOUTEROBJS=mld_zprecbld.o  mld_zprecset.o mld_zprecinit.o mld_zprecaply.o mld_zcprecset.o

COUTEROBJS=mld_cprecbld.o  mld_cprecset.o mld_cprecinit.o mld_cprecaply.o mld_ccprecset.o


OUTEROBJS=$(SOUTEROBJS) $(DOUTEROBJS) $(COUTEROBJS) $(ZOUTEROBJS)

F90OBJS=$(OUTEROBJS) $(INNEROBJS)

COBJS= mld_sslu_interface.o  \
       mld_dslu_interface.o mld_dumf_interface.o \
       mld_cslu_interface.o  \
       mld_zslu_interface.o mld_zumf_interface.o

OBJS=$(F90OBJS) $(COBJS) $(MPCOBJS)

LIBNAME=libmld_prec.a

lib:  $(OBJS) levd smoothd solvd
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)

levd:
	cd level && $(MAKE)
smoothd:
	cd smoother && $(MAKE)
solvd:
	cd solver && $(MAKE)

mpobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" F90COPT="$(F90COPT)")
	(make $(MPCOBJS) CC="$(MPCC)" CCOPT="$(CCOPT)")

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean: solvclean smoothclean levclean
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

levclean:
		cd level && $(MAKE) clean
smoothclean:
		cd smoother && $(MAKE) clean
solvclean:
		cd solver && $(MAKE) clean
