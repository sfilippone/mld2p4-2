include Make.inc

all: library 

library: libdir mlp kryl

libdir:
	(if test ! -d lib ; then mkdir lib; fi)
mlp:
	(cd mlprec; make lib)
kryl:
	(cd krylov; make symlink)
	(cd krylov; make lib)

install:
	($(INSTALL) -d $(INSTALL_DIR)/lib &&\
	   $(INSTALL_DATA) lib/*.a  $(INSTALL_DIR)/lib)
	($(INSTALL) -d $(INSTALL_DIR)/include && \
	   $(INSTALL_DATA) lib/*$(.mod) $(INSTALL_DIR)/include)
veryclean: 
	(cd mlprec; make veryclean)
	(cd krylov; make veryclean)
	(cd lib; /bin/rm -f *.a *$(.mod))

clean:
	(cd mlprec; make clean)
	(cd krylov; make clean)
