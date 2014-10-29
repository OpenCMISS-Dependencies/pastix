all:	src/_pastix.pxd \
	src/_murge.pxd	\
	src/_smurge.pxd	\
	src/_dmurge.pxd	\
	src/_cmurge.pxd	\
	src/_zmurge.pxd	\
	src/PASTIX.pxd	\
	src/_cscd_utils.pxd
	python setup.py install --prefix=${PASTIX_PREFIX}

test:
	-python -m doctest -v test_queue.py

clean:
	-rm -rf build
	-rm src/pypastix.c
