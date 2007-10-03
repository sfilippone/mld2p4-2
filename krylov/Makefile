include ../Make.inc
LIBDIR=../lib
PSBLIBDIR=$(PSBLASDIR)/lib
HERE=.
#
# Don't you dare touch the include order here!
# The method source file are symlinked, so as not to duplicate
# code, and the symlinks are recompiled using the local version
# of psb_prec_mod, so that they are tricked into using the
# MLD preconditioners while believing to be using the PSB ones.
#
FINCLUDES=$(FMFLAG). $(FMFLAG)$(LIBDIR) $(FMFLAG)$(PSBLIBDIR)
PSBKRYLDIR=$(PSBLASDIR)/krylov

METHDOBJS=psb_dcgstab.o psb_dcg.o psb_dcgs.o \
	psb_dbicg.o psb_dcgstabl.o psb_dgmresr.o\
	psb_zcgstab.o psb_zcgs.o psb_zgmresr.o

LIBMOD=psb_krylov_mod$(.mod)
MODOBJS=$(LIBMOD:$(.mod)=.o)
LOCAL_MODS=psb_prec_mod$(.mod)
LOCAL_OBJS=$(LOCAL_MODS:$(.mod)=.o)
OBJS=$(METHDOBJS) $(MODOBJS)
LIBNAME=libmld_krylov.a

lib: $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(LIBMOD) $(LIBDIR)

$(METHDOBJS): psb_prec_mod.o

symlink:
	(/bin/ln -fs $(PSBKRYLDIR)/*.f90 . )
	(/bin/ln -fs $(PSBKRYLDIR)/*.F90 . )

mld_krylov_mod.o: mld_krylov_mod.f90 $(LIBDIR)/mld_prec_mod$(.mod)
	$(F90) $(F90COPT) $(FINCLUDES) -c $<
psb_krylov_mod.o: psb_krylov_mod.f90 $(LIBDIR)/mld_prec_mod$(.mod)
	$(F90) $(F90COPT) $(FINCLUDES) -c $<


veryclean: clean
	(sh ./cleansymlinks *f90 *F90)
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS) $(LIBMOD) $(LOCAL_OBJS)
