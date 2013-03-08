include ../../../Make.inc
LIBDIR=../../../lib
PSBLIBDIR=$(PSBLASDIR)/lib
PSBINCDIR=$(PSBLASDIR)/include
HERE=../..

FINCLUDES=$(FMFLAG)../.. $(FMFLAG)$(LIBDIR) $(FMFLAG)$(PSBINCDIR) $(FMFLAG)$(PSBLIBDIR)


OBJS=mld_c_as_smoother_apply.o \
mld_c_as_smoother_apply_vect.o \
mld_c_as_smoother_bld.o \
mld_c_as_smoother_check.o \
mld_c_as_smoother_dmp.o \
mld_c_as_smoother_free.o \
mld_c_as_smoother_setc.o \
mld_c_as_smoother_seti.o \
mld_c_as_smoother_csetc.o \
mld_c_as_smoother_cseti.o \
mld_c_base_smoother_apply.o \
mld_c_base_smoother_apply_vect.o \
mld_c_base_smoother_bld.o \
mld_c_base_smoother_check.o \
mld_c_base_smoother_descr.o \
mld_c_base_smoother_dmp.o \
mld_c_base_smoother_free.o \
mld_c_base_smoother_setc.o \
mld_c_base_smoother_seti.o \
mld_c_base_smoother_setr.o \
mld_c_base_smoother_csetc.o \
mld_c_base_smoother_cseti.o \
mld_c_base_smoother_csetr.o \
mld_c_jac_smoother_apply.o \
mld_c_jac_smoother_apply_vect.o \
mld_c_jac_smoother_bld.o \
mld_d_as_smoother_apply.o \
mld_d_as_smoother_apply_vect.o \
mld_d_as_smoother_bld.o \
mld_d_as_smoother_check.o \
mld_d_as_smoother_dmp.o \
mld_d_as_smoother_free.o \
mld_d_as_smoother_setc.o \
mld_d_as_smoother_seti.o \
mld_d_as_smoother_csetc.o \
mld_d_as_smoother_cseti.o \
mld_d_base_smoother_apply.o \
mld_d_base_smoother_apply_vect.o \
mld_d_base_smoother_bld.o \
mld_d_base_smoother_check.o \
mld_d_base_smoother_descr.o \
mld_d_base_smoother_dmp.o \
mld_d_base_smoother_free.o \
mld_d_base_smoother_setc.o \
mld_d_base_smoother_seti.o \
mld_d_base_smoother_setr.o \
mld_d_base_smoother_csetc.o \
mld_d_base_smoother_cseti.o \
mld_d_base_smoother_csetr.o \
mld_d_jac_smoother_apply.o \
mld_d_jac_smoother_apply_vect.o \
mld_d_jac_smoother_bld.o \
mld_s_as_smoother_apply.o \
mld_s_as_smoother_apply_vect.o \
mld_s_as_smoother_bld.o \
mld_s_as_smoother_check.o \
mld_s_as_smoother_dmp.o \
mld_s_as_smoother_free.o \
mld_s_as_smoother_setc.o \
mld_s_as_smoother_seti.o \
mld_s_as_smoother_csetc.o \
mld_s_as_smoother_cseti.o \
mld_s_base_smoother_apply.o \
mld_s_base_smoother_apply_vect.o \
mld_s_base_smoother_bld.o \
mld_s_base_smoother_check.o \
mld_s_base_smoother_descr.o \
mld_s_base_smoother_dmp.o \
mld_s_base_smoother_free.o \
mld_s_base_smoother_setc.o \
mld_s_base_smoother_seti.o \
mld_s_base_smoother_setr.o \
mld_s_base_smoother_csetc.o \
mld_s_base_smoother_cseti.o \
mld_s_base_smoother_csetr.o \
mld_s_jac_smoother_apply.o \
mld_s_jac_smoother_apply_vect.o \
mld_s_jac_smoother_bld.o \
mld_z_as_smoother_apply.o \
mld_z_as_smoother_apply_vect.o \
mld_z_as_smoother_bld.o \
mld_z_as_smoother_check.o \
mld_z_as_smoother_dmp.o \
mld_z_as_smoother_free.o \
mld_z_as_smoother_setc.o \
mld_z_as_smoother_seti.o \
mld_z_as_smoother_csetc.o \
mld_z_as_smoother_cseti.o \
mld_z_base_smoother_apply.o \
mld_z_base_smoother_apply_vect.o \
mld_z_base_smoother_bld.o \
mld_z_base_smoother_check.o \
mld_z_base_smoother_descr.o \
mld_z_base_smoother_dmp.o \
mld_z_base_smoother_free.o \
mld_z_base_smoother_setc.o \
mld_z_base_smoother_seti.o \
mld_z_base_smoother_setr.o \
mld_z_base_smoother_csetc.o \
mld_z_base_smoother_cseti.o \
mld_z_base_smoother_csetr.o \
mld_z_jac_smoother_apply.o \
mld_z_jac_smoother_apply_vect.o \
mld_z_jac_smoother_bld.o 



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


