#----------------------------- tao makefile -----------------------------
#
#Read and modify options in ./src/Makefile
#
#Type  'make'  in this directory to compile tao 
#
#tao has been succesfully compiled with this Makefile in: 
# �iOS 11
# �IBM AIX Version 3.2 for IBM RISC 6000 workstations.
# �Linux for pentium processor. 
# �Hewlett Packard Envizex.
# �Sun Solaris OS5
#------------------------------------------------------------------------

include config.mk

all:
	(cd src; make)
	@echo; echo; echo Compilation done.
	@(echo "ADD ./tao/bin/ AND ./tao/script/ TO YOUR PATH.")
	@(echo "ADD  setenv tao_dir `pwd`  TO YOUR VARIABLES.")

clean_for_tar:
	(cd src; make clean)
	rm -f src/*.o lib/libreria.o lib/libreria.a 

tao: 
	(cd src; make tao)


vers: 	clean_for_tar
	rm -R -f tao tao_version
	mkdir tao tao/bin
	cp -R -L Makefile config.mk README demo doc include lib script src   tao
	rm -f tao/doc/first_compilation.txt #tao/lib/sistbanda* version_tmp/lib/surf_proc* version_tmp/lib/thin_sheet*
	tar -chf tao.tar tao
	gzip -f tao.tar
	echo "UPLOADING to github."
	touch tao/bin/touch_something #needed by git add
	mv tao tao_version
	(cd tao_version; git init; git remote add tao https://github.com/danigeos/tao-geo; git add .; git commit -a -mnewVersion; git push -u -f tao master)
