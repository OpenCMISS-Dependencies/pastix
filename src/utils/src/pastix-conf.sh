#!/bin/sh

LIBS="_LIBS_";
INC="_INC_";
CC="_CC_";
CXX="_CXX_"
CCOPT="_CCOPT_";
CXXOPT="_CXXOPT_";
CL="_CL_";
FC="_FC_";
FCOPT="_FCOPT_";
FL="_FL_";
OPTS="_OPTS_";
BLAS="_BLAS_";
VERSION="_VERSION_";
LIBS_MURGE=`echo $LIBS | sed -e 's/-lpastix/-lpastix_murge -lpastix/g'`

usage="usage : $0 [options] - Shows PaStiX libs, includes and compiler\n";
usage="$usage    options : \n";
usage="$usage        --libs               - prints librairies\n";
usage="$usage        --libs_murge         - prints librairies\n";
usage="$usage        --incs               - prints includes\n";
usage="$usage        --cc                 - prints C compiler\n";
usage="$usage        --ccopts             - prints C compiler options\n";
usage="$usage        --cxx                - prints C++ compiler\n";
usage="$usage        --cxxcopts           - prints C++ compiler options\n";
usage="$usage        --cl                 - prints C linker\n";
usage="$usage        --fc                 - prints fortran compiler\n";
usage="$usage        --fcopts             - prints fortran compiler options\n";
usage="$usage        --fl                 - prints fortran linker\n";
usage="$usage        --opts               - prints PaStiX compiling options\n";
usage="$usage        --vers               - prints PaStiX version\n";
usage="$usage        --blas               - prints blas choosen in config.in\n";

if [ $# = 0 ]
then
    echo "Librairies               : $LIBS" ;
    echo "Librairies with murge    : $LIBS_MURGE";
    echo "Incs                     : $INC" ;
    echo "C Compiler               : $CC" ;
    echo "C Compiler options       : $CCOPT" ;
    echo "C++ Compiler             : $CXX" ;
    echo "C++ Compiler options     : $CXXOPT" ;
    echo "Fortran Compiler         : $FC" ;
    echo "Fortran Compiler options : $FCOPT" ;
    echo "C Linker                 : $CL" ;
    echo "Fortran Linker           : $FL" ;
    echo "Options                  : $OPTS" ;
    echo "Version                  : $VERSION" ;
    echo "Blas                     : $BLAS" ;
elif [ $# = 1 ]
then
    case $1 in
        --libs)
            echo $LIBS;;
        --libs_murge)
            echo $LIBS_MURGE;;
        --incs)
            echo $INC;;
        --cc)
            echo $CC;;
        --ccopts)
            echo $CCOPT;;
        --cxx)
            echo $CXX;;
        --cxxopts)
            echo $CXXOPT;;
        --fc)
            echo $FC;;
        --fcopts)
            echo $FCOPT;;
        --cl)
            echo $CL;;
        --fl)
            echo $FL;;
        --opts)
            echo $OPTS;;
        --blas)
            echo $BLAS;;
        --vers)
            echo $VERSION;;

        *)
            echo -e $usage
    esac;
else
    echo -e $usage
fi;
