
from distutils.core import setup
from distutils.extension import Extension
import os
import fnmatch
import mpi4py

def touch(fname, times=None):
    with file(fname, 'a'):
        os.utime(fname, times)


def filter_filenames( filenames, match_list ):
    matched = []
    for name in filenames:
        for p in match_list:
            if fnmatch.fnmatch(name, p):
                matched.append( name )
                break # patterns
    return matched


def find_files(search_root, match_list):
    found = []
    for dirpath, dirnames, filenames in os.walk(search_root, followlinks = True):
        for d in dirnames[:]:
            if d.startswith("."):
                dirnames.remove(d)
        cur_dir_abs = os.path.abspath(dirpath)
        filenames_abs = []
        for name in filenames:
            if os.path.isabs(name):
                filenames_abs.append(name)
            else:
                filenames_abs.append(os.path.join(cur_dir_abs, name))
        found += filter_filenames(filenames_abs, match_list)
    return found



have_cython = False
try:
    from Cython.Distutils import build_ext as _build_ext
    have_cython = True
except ImportError:
    from distutils.command.build_ext import build_ext as _build_ext
    print("*"*80)
    print()
    print(" You need to generate C source files with Cython!!")
    print(" Download and install Cython <http://www.cython.org>")
    print()
    print("*"*80)

class build_ext(_build_ext):

    user_options = _build_ext.user_options + [ ('include_dirs=', None, 'include directories'),
                                               ('library_dirs=', None, 'library directories'),
                                               ('libraries=', None, 'libraries')
                                             ]

    def finalize_options(self):
        self.ensure_string_list('libraries')
        _build_ext.finalize_options(self)


import numpy

if have_cython:
    touch('src/pypastix.pyx')
    ext_pypastix  = Extension('pypastix', ['src/pypastix.pyx'],
                              extra_compile_args = ["-O3", "-Wall"],
                              include_dirs = [numpy.get_include(),
                                              mpi4py.get_include()])
else:
    touch('src/pypastix.c')
    ext_pypastix  = Extension('pypastix', ['src/pypastix.c'],
                              extra_compile_args = ["-O3", "-Wall"],
                              include_dirs = [numpy.get_include(),
                                              mpi4py.get_include()],
                          )


setup(
    name = 'pypastix',
    version='0.0.0',
    packages= [],
    author='Xavier LACOSTE',
    author_email='xl64100@gmail.com',
    description='Cython wrapper around PaStiX',
    keywords = ["pastix", "cython", "sparse matrix", "hpc"],
    license = 'bsd',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [ext_pypastix]
)
