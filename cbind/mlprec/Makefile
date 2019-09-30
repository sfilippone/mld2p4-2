include ../../Make.inc
LIBDIR=../../lib
INCDIR=../../include
MODDIR=../../modules
HERE=.

CINCLUDES=-I. -I$(LIBDIR) -I$(PSBLAS_INCDIR)
FINCLUDES=$(FMFLAG)$(HERE) $(FMFLAG)$(INCDIR) $(FMFLAG)$(MODDIR) $(PSBLAS_INCLUDES)


OBJS=mld_prec_cbind_mod.o mld_dprec_cbind_mod.o  mld_c_dprec.o
CMOD=mld_cbind.h mld_c_dprec.h mld_const.h


LIBMOD=mld_prec_cbind_mod$(.mod) 
LOCAL_MODS=$(LIBMOD)
#LIBNAME=$(CPRECLIBNAME)	


lib:  $(OBJS) $(CMOD)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(LIBMOD) $(CMOD) $(INCDIR)

mld_prec_cbind_mod.o: mld_dprec_cbind_mod.o 
#mld_prec_cbind_mod.o: psb_prec_cbind_mod.o
veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
