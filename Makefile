include Make.inc


all:  library 

library: libdir mlp 

libdir:
	(if test ! -d lib ; then mkdir lib; fi)
	(if test ! -d include ; then mkdir include; fi)
	($(INSTALL_DATA) Make.inc  include/Make.inc.mld2p4)
        

mlp:
	cd mlprec && $(MAKE) all

install: all
	(./mkdir.sh  $(INSTALL_DIR) &&\
	   $(INSTALL_DATA) Make.inc  $(INSTALL_DIR)/Make.inc.MLD2P4)
	(./mkdir.sh $(INSTALL_LIBDIR) &&\
	   $(INSTALL_DATA) lib/*.a  $(INSTALL_LIBDIR))
	(./mkdir.sh $(INSTALL_INCLUDEDIR) && \
	   $(INSTALL_DATA) include/*$(.mod) $(INSTALL_INCLUDEDIR))
	(./mkdir.sh  $(INSTALL_INCLUDEDIR) &&\
	   $(INSTALL_DATA) Make.inc  $(INSTALL_INCLUDEDIR)/Make.inc.mld2p4)
	(./mkdir.sh $(INSTALL_INCLUDEDIR) && \
	   $(INSTALL_DATA) include/*.h $(INSTALL_INCLUDEDIR))
	(./mkdir.sh  $(INSTALL_DOCSDIR) && \
	   /bin/cp -fr docs/*pdf docs/html $(INSTALL_DOCSDIR))
	(./mkdir.sh  $(INSTALL_DOCSDIR) && \
	   $(INSTALL_DATA) README LICENSE $(INSTALL_DOCSDIR))
cleanlib:
	(cd lib; /bin/rm -f *.a *$(.mod) *$(.fh))
	(cd include; /bin/rm -f *.a *$(.mod) *$(.fh))

veryclean: cleanlib
	(cd mlprec; make veryclean)
	(cd examples/fileread; make clean)
	(cd examples/pdegen; make clean)
	(cd tests/fileread; make clean)
	(cd tests/pdegen; make clean)

check: all
	make check -C tests/pdegen

clean:
	(cd mlprec; make clean)
