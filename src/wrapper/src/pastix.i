%module pastix
%{
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifndef FORCE_NOMPI
#include <mpi.h>
#endif
#include <complex.h>
#include "pastix.h"
%}
%include "cpointer.i"
%include "carrays.i"
%include "mpi4py.i"

%mpi4py_typemap(Comm, MPI_Comm);

%typemap(in) int32_t {
   $1 = PyInt_AsLong($input);
 }
%typemap(in) int64_t {
   $1 = PyInt_AsLong($input);
 }

%array_functions(int, intArray);
%array_functions(pastix_int_t, pIntArray);
%array_functions(double, doubleArray);
%array_functions(float, simpleArray);
%array_functions(pastix_float_t, pFloatArray);

%apply long* INOUT{pastix_data_t **pastix_data }; // tell SWIG its an inout
#%include "pastix.h"
