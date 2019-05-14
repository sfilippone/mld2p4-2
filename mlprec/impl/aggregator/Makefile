include ../../../Make.inc
LIBDIR=../../../lib
INCDIR=../../../include
MODDIR=../../../modules 
HERE=../..

FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(PSBLAS_INCLUDES)

#CINCLUDES= -I${SUPERLU_INCDIR} -I${HSL_INCDIR}  -I${SPRAL_INCDIR} -I/home/users/pasqua/Ambra/BootCMatch/include -lBCM -L/home/users/pasqua/Ambra/BootCMatch/lib  -lm 

OBJS= \
mld_s_dec_aggregator_mat_asb.o \
mld_s_dec_aggregator_mat_bld.o \
mld_s_dec_aggregator_tprol.o \
mld_s_symdec_aggregator_tprol.o  \
mld_s_map_to_tprol.o mld_s_soc1_map_bld.o mld_s_soc2_map_bld.o\
mld_saggrmat_biz_bld.o mld_saggrmat_minnrg_bld.o\
mld_saggrmat_nosmth_bld.o mld_saggrmat_smth_bld.o \
mld_d_dec_aggregator_mat_asb.o \
mld_d_dec_aggregator_mat_bld.o \
mld_d_dec_aggregator_tprol.o \
mld_d_symdec_aggregator_tprol.o  \
mld_d_map_to_tprol.o mld_d_soc1_map_bld.o mld_d_soc2_map_bld.o\
mld_daggrmat_biz_bld.o mld_daggrmat_minnrg_bld.o\
mld_daggrmat_nosmth_bld.o mld_daggrmat_smth_bld.o \
mld_c_dec_aggregator_mat_asb.o \
mld_c_dec_aggregator_mat_bld.o \
mld_c_dec_aggregator_tprol.o \
mld_c_symdec_aggregator_tprol.o  \
mld_c_map_to_tprol.o mld_c_soc1_map_bld.o mld_c_soc2_map_bld.o\
mld_caggrmat_biz_bld.o mld_caggrmat_minnrg_bld.o\
mld_caggrmat_nosmth_bld.o mld_caggrmat_smth_bld.o \
mld_z_dec_aggregator_mat_asb.o \
mld_z_dec_aggregator_mat_bld.o \
mld_z_dec_aggregator_tprol.o \
mld_z_symdec_aggregator_tprol.o  \
mld_z_map_to_tprol.o mld_z_soc1_map_bld.o mld_z_soc2_map_bld.o\
mld_zaggrmat_biz_bld.o mld_zaggrmat_minnrg_bld.o\
mld_zaggrmat_nosmth_bld.o mld_zaggrmat_smth_bld.o 

#mld_s_hybrid_aggregator_tprol.o \
#mld_d_hybrid_aggregator_tprol.o \
#mld_c_hybrid_aggregator_tprol.o \
#mld_z_hybrid_aggregator_tprol.o \

#bootCMatch_interface.o mld_d_bcmatch_aggregator_tprol.o\
#mld_d_bcmatch_map_to_tprol.o mld_d_bcmatch_aggregator_mat_bld.o \


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
