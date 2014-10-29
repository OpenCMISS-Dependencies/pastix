Topic: To use PaStiX in Python

       You can customize the instalation PATH by setting PYTHON_PREFIX in
       config.in file, the directory $PYTHON_PREFIX/pythonX.Y/site-packages
       must be in your PYTHONPATH at runtime.

       To use PaStiX with Python you have to install Cython, then simply run:
       >  make pypastix
       And test it with
       > python ./example/src/pypastix/pastix.py

