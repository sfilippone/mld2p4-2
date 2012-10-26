include ../../Make.inc
LIBDIR=../../lib
PSBLIBDIR=$(PSBLASDIR)/lib
PSBINCDIR=$(PSBLASDIR)/include
HERE=..

FINCLUDES=$(FMFLAG).. $(FMFLAG)$(LIBDIR) $(FMFLAG)$(PSBINCDIR) $(FMFLAG)$(PSBLIBDIR)


DMPFOBJS=mld_daggrmat_nosmth_asb.o mld_daggrmat_smth_asb.o mld_daggrmat_minnrg_asb.o mld_daggrmat_biz_asb.o

SMPFOBJS=mld_saggrmat_nosmth_asb.o mld_saggrmat_smth_asb.o mld_saggrmat_minnrg_asb.o mld_saggrmat_biz_asb.o

ZMPFOBJS=mld_zaggrmat_nosmth_asb.o mld_zaggrmat_smth_asb.o mld_zaggrmat_minnrg_asb.o mld_zaggrmat_biz_asb.o

CMPFOBJS=mld_caggrmat_nosmth_asb.o mld_caggrmat_smth_asb.o mld_caggrmat_minnrg_asb.o mld_caggrmat_biz_asb.o


MPFOBJS=$(SMPFOBJS) $(DMPFOBJS) $(CMPFOBJS) $(ZMPFOBJS) 

MPCOBJS=mld_sslud_interface.o mld_dslud_interface.o mld_cslud_interface.o mld_zslud_interface.o


DINNEROBJS= mld_dcoarse_bld.o  mld_dmlprec_bld.o	  \
	mld_dilu0_fact.o   mld_diluk_fact.o  mld_dilut_fact.o  mld_daggrmap_bld.o \
	mld_d_dec_map_bld.o mld_dmlprec_aply.o mld_daggrmat_asb.o   \
	$(DMPFOBJS)

# mld_d_base_solver_impl.o mld_d_base_smoother_impl.o mld_d_onelev_impl.o\
	mld_d_as_smoother_impl.o mld_d_jac_smoother_impl.o \
	mld_d_diag_solver_impl.o mld_d_id_solver_impl.o mld_d_ilu_solver_impl.o

SINNEROBJS= mld_scoarse_bld.o  mld_smlprec_bld.o	  \
	mld_silu0_fact.o   mld_siluk_fact.o  mld_silut_fact.o  mld_saggrmap_bld.o \
	mld_s_dec_map_bld.o mld_smlprec_aply.o mld_saggrmat_asb.o   \
	$(SMPFOBJS)

#  mld_s_base_solver_impl.o mld_s_base_smoother_impl.o mld_s_onelev_impl.o\
	mld_s_as_smoother_impl.o mld_s_jac_smoother_impl.o \
	mld_s_diag_solver_impl.o mld_s_id_solver_impl.o mld_s_ilu_solver_impl.o

ZINNEROBJS= mld_zcoarse_bld.o  mld_zmlprec_bld.o	  \
	mld_zilu0_fact.o   mld_ziluk_fact.o  mld_zilut_fact.o  mld_zaggrmap_bld.o \
	mld_z_dec_map_bld.o mld_zmlprec_aply.o mld_zaggrmat_asb.o   \
	$(ZMPFOBJS) 

# mld_z_base_solver_impl.o mld_z_base_smoother_impl.o mld_z_onelev_impl.o\
	mld_z_as_smoother_impl.o mld_z_jac_smoother_impl.o \
	mld_z_diag_solver_impl.o mld_z_id_solver_impl.o mld_z_ilu_solver_impl.o

CINNEROBJS= mld_ccoarse_bld.o  mld_cmlprec_bld.o	  \
	mld_cilu0_fact.o   mld_ciluk_fact.o  mld_cilut_fact.o  mld_caggrmap_bld.o \
	mld_c_dec_map_bld.o mld_cmlprec_aply.o mld_caggrmat_asb.o   \
	$(CMPFOBJS) 

# mld_c_base_solver_impl.o mld_c_base_smoother_impl.o mld_c_onelev_impl.o\
	mld_c_as_smoother_impl.o mld_c_jac_smoother_impl.o \
	mld_c_diag_solver_impl.o mld_c_id_solver_impl.o mld_c_ilu_solver_impl.o

INNEROBJS= $(SINNEROBJS) $(DINNEROBJS) $(CINNEROBJS) $(ZINNEROBJS)


DOUTEROBJS=mld_dprecbld.o  mld_dprecset.o mld_dprecinit.o mld_dprecaply.o

SOUTEROBJS=mld_sprecbld.o  mld_sprecset.o mld_sprecinit.o mld_sprecaply.o

ZOUTEROBJS=mld_zprecbld.o  mld_zprecset.o mld_zprecinit.o mld_zprecaply.o

COUTEROBJS=mld_cprecbld.o  mld_cprecset.o mld_cprecinit.o mld_cprecaply.o


OUTEROBJS=$(SOUTEROBJS) $(DOUTEROBJS) $(COUTEROBJS) $(ZOUTEROBJS)

F90OBJS=$(OUTEROBJS) $(INNEROBJS)

COBJS= mld_sslu_interface.o mld_sumf_interface.o \
       mld_dslu_interface.o mld_dumf_interface.o \
       mld_cslu_interface.o mld_cumf_interface.o \
       mld_zslu_interface.o mld_zumf_interface.o

OBJS=$(F90OBJS) $(COBJS) $(MPCOBJS)

LIBNAME=libmld_prec.a

lib:  $(OBJS) smoothd solvd
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)

smoothd:
	cd SMOOTHER && $(MAKE)
solvd:
	cd SOLVER && $(MAKE)

mpobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" F90COPT="$(F90COPT)")
	(make $(MPCOBJS) CC="$(MPCC)" CCOPT="$(CCOPT)")

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean: solvclean smoothclean
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

smoothclean:
		cd SMOOTHER && $(MAKE) clean
solvclean:
		cd SOLVER && $(MAKE) clean
