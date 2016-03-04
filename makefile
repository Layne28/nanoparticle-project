np:
	make -f makefile.np

etch:
	make -f makefile.etch

irrev:

	make -f makefile.irrev

glauber:

	make -f makefile.glauber

np_cubic:

	make -f makefile.np_cubic

kmc:

	make -f makefile.kmc

etch_cubic:

	make -f makefile.etch_cubic

.PHONY:analyze
analyze:

	make -f makefile.analyze

clean:

	rm *.o *.out
