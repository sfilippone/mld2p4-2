include ../Make.inc
LIBDIR=../lib
PSBLIBDIR=$(PSBLASDIR)/lib
PSBINCDIR=$(PSBLASDIR)/include
HERE=.

FINCLUDES=$(FMFLAG). $(FMFLAG)$(LIBDIR) $(FMFLAG)$(PSBINCDIR) $(FMFLAG)$(PSBLIBDIR)


SMODOBJS=mld_s_prec_type.o mld_s_prec_mod.o mld_s_move_alloc_mod.o mld_s_ilu_fact_mod.o \
	mld_s_inner_mod.o mld_s_ilu_solver.o mld_s_diag_solver.o mld_s_jac_smoother.o mld_s_as_smoother.o \
	mld_s_id_solver.o mld_s_slu_solver.o

DMODOBJS=mld_d_prec_type.o mld_d_prec_mod.o mld_d_move_alloc_mod.o mld_d_ilu_fact_mod.o \
	mld_d_inner_mod.o mld_d_ilu_solver.o mld_d_diag_solver.o mld_d_jac_smoother.o mld_d_as_smoother.o \
	mld_d_umf_solver.o mld_d_slu_solver.o mld_d_sludist_solver.o mld_d_id_solver.o 

CMODOBJS=mld_c_prec_type.o mld_c_prec_mod.o mld_c_move_alloc_mod.o mld_c_ilu_fact_mod.o \
	mld_c_inner_mod.o mld_c_ilu_solver.o mld_c_diag_solver.o mld_c_jac_smoother.o mld_c_as_smoother.o \
	mld_c_id_solver.o mld_c_slu_solver.o

ZMODOBJS=mld_z_prec_type.o mld_z_prec_mod.o mld_z_move_alloc_mod.o mld_z_ilu_fact_mod.o \
	mld_z_inner_mod.o mld_z_ilu_solver.o mld_z_diag_solver.o mld_z_jac_smoother.o mld_z_as_smoother.o \
	mld_z_id_solver.o mld_z_umf_solver.o mld_z_slu_solver.o 

MODOBJS=mld_base_prec_type.o mld_prec_type.o  mld_prec_mod.o \
	$(SMODOBJS) $(DMODOBJS) $(CMODOBJS) $(ZMODOBJS) 

SMPFOBJS=mld_saggrmat_nosmth_asb.o mld_saggrmat_smth_asb.o 
DMPFOBJS=mld_daggrmat_nosmth_asb.o mld_daggrmat_smth_asb.o mld_daggrmat_minnrg_asb.o 
CMPFOBJS=mld_caggrmat_nosmth_asb.o mld_caggrmat_smth_asb.o 
ZMPFOBJS=mld_zaggrmat_nosmth_asb.o mld_zaggrmat_smth_asb.o 
MPFOBJS=$(SMPFOBJS) $(DMPFOBJS) $(CMPFOBJS) $(ZMPFOBJS) 

MPCOBJS=mld_sslud_interface.o mld_dslud_interface.o mld_cslud_interface.o mld_zslud_interface.o

SINNEROBJS= mld_scoarse_bld.o  mld_smlprec_bld.o     \
	mld_silu0_fact.o   mld_siluk_fact.o  mld_silut_fact.o  mld_saggrmap_bld.o \
	mld_s_dec_map_bld.o mld_smlprec_aply.o mld_saggrmat_asb.o   \
	$(SMPFOBJS) 

DINNEROBJS= mld_dcoarse_bld.o  mld_dmlprec_bld.o	  \
	mld_dilu0_fact.o   mld_diluk_fact.o  mld_dilut_fact.o  mld_daggrmap_bld.o \
	mld_d_dec_map_bld.o mld_dmlprec_aply.o mld_daggrmat_asb.o   \
	$(DMPFOBJS) 

CINNEROBJS= mld_ccoarse_bld.o  mld_cmlprec_bld.o \
	mld_cilu0_fact.o   mld_ciluk_fact.o  mld_cilut_fact.o  mld_caggrmap_bld.o \
	mld_c_dec_map_bld.o mld_cmlprec_aply.o mld_caggrmat_asb.o   \
	$(CMPFOBJS) 

ZINNEROBJS= mld_zcoarse_bld.o  mld_zmlprec_bld.o \
	mld_zilu0_fact.o   mld_ziluk_fact.o  mld_zilut_fact.o  mld_zaggrmap_bld.o \
	mld_z_dec_map_bld.o mld_zmlprec_aply.o mld_zaggrmat_asb.o   \
	$(ZMPFOBJS) 

INNEROBJS= $(SINNEROBJS) $(DINNEROBJS) $(CINNEROBJS) $(ZINNEROBJS)

SOUTEROBJS=mld_sprecbld.o  mld_sprecset.o mld_sprecinit.o mld_sprecaply.o 
DOUTEROBJS=mld_dprecbld.o  mld_dprecset.o mld_dprecinit.o mld_dprecaply.o
COUTEROBJS=mld_cprecbld.o  mld_cprecset.o mld_cprecinit.o mld_cprecaply.o 
ZOUTEROBJS=mld_zprecbld.o  mld_zprecset.o mld_zprecinit.o mld_zprecaply.o 

OUTEROBJS=$(SOUTEROBJS) $(DOUTEROBJS) $(COUTEROBJS) $(ZOUTEROBJS)

F90OBJS=$(OUTEROBJS) $(INNEROBJS)

COBJS= mld_sslu_interface.o mld_sumf_interface.o \
       mld_dslu_interface.o mld_dumf_interface.o \
       mld_cslu_interface.o mld_cumf_interface.o \
       mld_zslu_interface.o mld_zumf_interface.o

OBJS=$(F90OBJS) $(MODOBJS) $(COBJS) $(MPCOBJS)

LIBMOD=mld_prec_mod$(.mod) 
LOCAL_MODS=$(MODOBJS:.o=$(.mod))
LIBNAME=libmld_prec.a

lib: $(LIBDIR)/$(LIBNAME)

$(LIBNAME): $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)


# flea: if libdir misses some .mod file, it won't be copied unless the .a file is missing too  
$(LIBDIR)/$(LIBNAME): $(LIBNAME) 
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(LIBMOD) $(LOCAL_MODS) mld_const.h $(LIBDIR)

mld_s_prec_type.o mld_d_prec_type.o mld_c_prec_type.o mld_z_prec_type.o	: mld_base_prec_type.o
mld_prec_type.o: mld_s_prec_type.o mld_d_prec_type.o mld_c_prec_type.o mld_z_prec_type.o	
mld_prec_mod.o: mld_prec_type.o mld_s_prec_mod.o mld_d_prec_mod.o mld_c_prec_mod.o mld_z_prec_mod.o

$(MODOBJS): $(PSBINCDIR)/psb_base_mod$(.mod)

$(SINNEROBJS) $(SOUTEROBJS): $(SMODOBJS)
$(DINNEROBJS) $(DOUTEROBJS): $(DMODOBJS)
$(CINNEROBJS) $(COUTEROBJS): $(CMODOBJS)
$(ZINNEROBJS) $(ZOUTEROBJS): $(ZMODOBJS)

mld_s_inner_mod.o: mld_s_move_alloc_mod.o mld_s_prec_type.o
mld_d_inner_mod.o: mld_d_move_alloc_mod.o mld_d_prec_type.o
mld_c_inner_mod.o: mld_c_move_alloc_mod.o mld_c_prec_type.o
mld_z_inner_mod.o: mld_z_move_alloc_mod.o mld_z_prec_type.o

mld_s_move_alloc_mod.o: mld_s_prec_type.o
mld_d_move_alloc_mod.o: mld_d_prec_type.o
mld_c_move_alloc_mod.o: mld_c_prec_type.o
mld_z_move_alloc_mod.o: mld_z_prec_type.o

mld_s_prec_mod.o: mld_s_move_alloc_mod.o 
mld_d_prec_mod.o: mld_d_move_alloc_mod.o 
mld_c_prec_mod.o: mld_c_move_alloc_mod.o 
mld_z_prec_mod.o: mld_z_move_alloc_mod.o 



mld_d_sludist_solver.o mld_d_slu_solver.o mld_d_umf_solver.o mld_d_diag_solver.o mld_d_ilu_solver.o: mld_d_prec_type.o 
mld_d_ilu_fact_mod.o: mld_base_prec_type.o	
mld_d_ilu_solver.o mld_d_iluk_fact.o: mld_d_ilu_fact_mod.o
mld_d_as_smoother.o mld_d_jac_smoother.o: mld_d_prec_type.o 
mld_d_jac_smoother.o: mld_d_diag_solver.o
mld_dprecinit.o mld_dprecset.o: mld_d_diag_solver.o mld_d_ilu_solver.o \
		mld_d_umf_solver.o mld_d_as_smoother.o mld_d_jac_smoother.o \
	        mld_d_id_solver.o mld_d_slu_solver.o mld_d_sludist_solver.o

mld_z_umf_solver.o mld_z_diag_solver.o mld_z_ilu_solver.o: mld_z_prec_type.o 
mld_z_ilu_fact_mod.o: mld_base_prec_type.o	
mld_z_ilu_solver.o mld_z_iluk_fact.o: mld_z_ilu_fact_mod.o
mld_z_as_smoother.o mld_z_jac_smoother.o: mld_z_prec_type.o 
mld_z_jac_smoother.o: mld_z_diag_solver.o
mld_zprecinit.o mld_zprecset.o: mld_z_diag_solver.o mld_z_ilu_solver.o \
		mld_z_umf_solver.o mld_z_as_smoother.o mld_z_jac_smoother.o

mld_s_diag_solver.o mld_s_ilu_solver.o: mld_s_prec_type.o 
mld_s_ilu_fact_mod.o: mld_base_prec_type.o	
mld_s_ilu_solver.o mld_s_iluk_fact.o: mld_s_ilu_fact_mod.o
mld_s_as_smoother.o mld_s_jac_smoother.o: mld_s_prec_type.o 
mld_s_jac_smoother.o: mld_s_diag_solver.o
mld_sprecinit.o mld_sprecset.o: mld_s_diag_solver.o mld_s_ilu_solver.o \
		mld_s_as_smoother.o mld_s_jac_smoother.o

mld_c_diag_solver.o mld_c_ilu_solver.o: mld_c_prec_type.o 
mld_c_ilu_fact_mod.o: mld_base_prec_type.o	
mld_c_ilu_solver.o mld_c_iluk_fact.o: mld_c_ilu_fact_mod.o
mld_c_as_smoother.o mld_c_jac_smoother.o: mld_c_prec_type.o 
mld_c_jac_smoother.o: mld_c_diag_solver.o
mld_cprecinit.o mld_cprecset.o: mld_c_diag_solver.o mld_c_ilu_solver.o \
		mld_c_as_smoother.o mld_c_jac_smoother.o

 
mpobjs: $(MODOBJS)
	(make $(MPFOBJS) F90="$(MPF90)" F90COPT="$(F90COPT)")
	(make $(MPCOBJS) CC="$(MPCC)" CCOPT="$(CCOPT)")

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

