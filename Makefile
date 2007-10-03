include Make.inc

library: mlp kryl

mlp:
	(cd mlprec; make lib)
kryl:
	(cd krylov; make symlink)
	(cd krylov; make lib)
veryclean: 
	(cd mlprec; make veryclean)
	(cd krylov; make veryclean)
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

clean:
	(cd mlprec; make clean)
	(cd krylov; make clean)
	/bin/rm -f $(OBJS) $(LOCAL_MODS)
