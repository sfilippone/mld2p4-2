include ../../Make.inc
LIBDIR=../../lib
INCDIR=../../include
MODDIR=../../modules 
HERE=..

FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(PSBLAS_INCLUDES)



DMPFOBJS=

SMPFOBJS=

ZMPFOBJS=

CMPFOBJS=


MPFOBJS=$(SMPFOBJS) $(DMPFOBJS) $(CMPFOBJS) $(ZMPFOBJS) 

MPCOBJS=mld_dslud_interface.o mld_zslud_interface.o


DINNEROBJS= mld_dmlprec_bld.o mld_dfile_prec_descr.o \
	mld_d_smoothers_bld.o mld_d_hierarchy_bld.o \
	mld_dmlprec_aply.o   \
	$(DMPFOBJS) mld_d_extprol_bld.o

SINNEROBJS=   mld_smlprec_bld.o mld_sfile_prec_descr.o \
	mld_s_smoothers_bld.o mld_s_hierarchy_bld.o \
	mld_smlprec_aply.o   \
	$(SMPFOBJS) mld_s_extprol_bld.o 

ZINNEROBJS=   mld_zmlprec_bld.o mld_zfile_prec_descr.o \
	mld_z_smoothers_bld.o mld_z_hierarchy_bld.o \
	mld_zmlprec_aply.o   \
	$(ZMPFOBJS)  mld_z_extprol_bld.o

CINNEROBJS=   mld_cmlprec_bld.o mld_cfile_prec_descr.o \
	mld_c_smoothers_bld.o mld_c_hierarchy_bld.o \
	mld_cmlprec_aply.o   \
	$(CMPFOBJS)  mld_c_extprol_bld.o

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

lib:  $(OBJS) aggrd levd smoothd solvd
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)

aggrd:
	cd aggregator && $(MAKE)
levd:
	cd level && $(MAKE)
smoothd:
	cd smoother && $(MAKE)
solvd:
	cd solver && $(MAKE)

mpobjs: 
	(make $(MPFOBJS) FC="$(MPFC)" FCOPT="$(FCOPT)")
	(make $(MPCOBJS) CC="$(MPCC)" CCOPT="$(CCOPT)")

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean: solvclean smoothclean levclean aggrclean
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

aggrclean:
		cd aggregator && $(MAKE) clean
levclean:
		cd level && $(MAKE) clean
smoothclean:
		cd smoother && $(MAKE) clean
solvclean:
		cd solver && $(MAKE) clean
