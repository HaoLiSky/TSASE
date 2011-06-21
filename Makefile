FCOMP=gnu95

all:
	python -c 'import compileall; compileall.compile_dir(".",force=True)' #quiet=True
	find calc -name 'Makefile' | xargs perl -pi -e 's/fcompiler=[\w]+/fcompiler=$(FCOMP)/g'
	cd calc && $(MAKE)

.PHONY: clean
clean:
	cd calc && $(MAKE) clean
	find . -name "*.pyc" -exec rm {} \;
