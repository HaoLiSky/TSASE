FCOMP=gnu95

all:
	python -c 'import compileall; compileall.compile_dir(".",force=True)' #quiet=True
	find calculators -name 'Makefile' | xargs perl -pi -e 's/fcompiler=[\w]+/fcompiler=$(FCOMP)/g'
	cd calculators && $(MAKE)

.PHONY: clean
clean:
	cd calculators && $(MAKE) clean
	find . -name "*.pyc" -exec rm {} \;
