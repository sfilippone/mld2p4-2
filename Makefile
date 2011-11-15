include Make.inc

all: library 

library: libdir mlp 

libdir:
	(if test ! -d lib ; then mkdir lib; fi)
mlp:
	cd mlprec && $(MAKE) lib

install:
	(./mkdir.sh  $(INSTALL_DIR) &&\
	   $(INSTALL_DATA) Make.inc  $(INSTALL_DIR)/Make.inc.MLD2P4)
	(./mkdir.sh $(INSTALL_LIBDIR) &&\
	   $(INSTALL_DATA) lib/*.a  $(INSTALL_LIBDIR))
	(./mkdir.sh $(INSTALL_INCLUDEDIR) && \
	   $(INSTALL_DATA) lib/*$(.mod) $(INSTALL_INCLUDEDIR))
	(./mkdir.sh $(INSTALL_INCLUDEDIR) && \
	   $(INSTALL_DATA) lib/*.h $(INSTALL_INCLUDEDIR))
	(./mkdir.sh  $(INSTALL_DOCSDIR) && \
	   /bin/cp -fr docs/*pdf docs/html $(INSTALL_DOCSDIR))
veryclean: 
	(cd mlprec; make veryclean)
	(cd lib; /bin/rm -f *.a *$(.mod))
	(cd examples/fileread; make clean)
	(cd examples/pdegen; make clean)
	(cd tests/fileread; make clean)
	(cd tests/pdegen; make clean)

check: all
	make check -C tests/pdegen

clean:
	(cd mlprec; make clean)
