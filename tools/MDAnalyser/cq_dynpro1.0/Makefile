sinclude make.sys

default :
	@echo '  make dynpro'

dynpro : bindir libqemod flibs 
	if test -d src ; then \
	( cd src ; $(MAKE) TLDEPS= all || exit 1) ; fi

libqemod :
	( cd modules ; $(MAKE) TLDEPS= all || exit 1 )
flibs :
	( cd flib ; $(MAKE) TLDEPS= $(FLIB_TARGETS) || exit 1 )

bindir :

#########################################################
# plugins
#########################################################

links : bindir
	( cd bin/ ; \
	rm -f *.x ; \
	for exe in ../*/*/*.x ../*/bin/* ; do \
	    if test ! -L $$exe ; then ln -fs $$exe . ; fi \
	done \
	)

bindir :
	test -d bin || mkdir bin

#########################################################
# Other targets: clean up
#########################################################

# remove object files and executables
clean :
	touch make.sys 
	for dir in \
		src modules flib \
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		$(MAKE) TLDEPS= clean ) \
	    fi \
	done
	- /bin/rm -rf bin/*.x tmp 
	- /bin/rm -rf src/*~ 
#	- cd examples/example01; /bin/rm -rf CRASH out.* pch.* sph.* UpdatedAtoms.* ; cd -

tar :
	@if test -f dynpro.tar.gz ; then /bin/rm dynpro.tar.gz ; fi
	# do not include unneeded stuff 
	find ./ -type f | grep -v -e /.svn/ -e'/\.' -e'\.o$$' \
             -e'\.mod$$' -e'\.a$$' -e'\.d$$' -e'\.i$$' -e'\.F90$$' -e'\.x$$' \
	     -e'~$$' -e'\./GUI' -e '\./tempdir' | xargs tar rvf dynpro.tar
	gzip dynpro.tar

