TOP=../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)/lib
INCDIR=$(TOP)/include
MODDIR=$(TOP)/modules
HERE=.
DEST=../

CINCLUDES=-I. -I$(LIBDIR) -I$(PSBLAS_INCDIR)
FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(INCDIR) $(FMFLAG)$(MODDIR) $(PSBLAS_INCLUDES)


OBJS=mld_prec_cbind_mod.o mld_dprec_cbind_mod.o  mld_c_dprec.o mld_zprec_cbind_mod.o  mld_c_zprec.o
CMOD=mld_cbind.h mld_c_dprec.h  mld_c_zprec.h mld_const.h


LIBMOD=mld_prec_cbind_mod$(.mod) 
LOCAL_MODS=$(LIBMOD)
#LIBNAME=$(CPRECLIBNAME)	


lib:  $(OBJS) $(CMOD)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(DEST)
	/bin/cp -p $(LIBMOD) $(CMOD) $(DEST)

mld_prec_cbind_mod.o: mld_dprec_cbind_mod.o  mld_zprec_cbind_mod.o 
#mld_prec_cbind_mod.o: psb_prec_cbind_mod.o
veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
