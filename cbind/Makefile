include ../Make.inc

HERE=.
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules/
LIBNAME=$(CBINDLIBNAME)
LIBNAME=libmld_cbind.a	

lib: mlprecd
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(CPUPDFLAG) *.h $(INCDIR)
	/bin/cp -p $(CPUPDFLAG) *$(.mod) $(MODDIR)	


mlprecd:
	cd mlprec && $(MAKE) lib LIBNAME=$(LIBNAME) 


clean: 
	cd mlprec && $(MAKE) clean


veryclean: clean
	cd test/pargen && $(MAKE) clean
	/bin/rm -f $(HERE)/$(LIBNAME) $(LIBMOD) *$(.mod) *.h
