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
	(./mkdir.sh $(INSTALL_LIBDIR) &&\
	   $(INSTALL_DATA) lib/*.a  $(INSTALL_LIBDIR))
	(./mkdir.sh $(INSTALL_INCLUDEDIR) && \
	   $(INSTALL_DATA) lib/*$(.mod) $(INSTALL_INCLUDEDIR))
	(./mkdir.sh  $(INSTALL_DOCSDIR) && \
	   /bin/cp -fr docs/*pdf docs/html $(INSTALL_DOCSDIR))
veryclean: 
	(cd mlprec; make veryclean)
	(cd krylov; make veryclean)
	(cd lib; /bin/rm -f *.a *$(.mod))
	(cd examples/fileread; make clean)
	(cd examples/pdegen; make clean)
	(cd tests/fileread; make clean)
	(cd tests/pdegen; make clean)

clean:
	(cd mlprec; make clean)
	(cd krylov; make clean)
