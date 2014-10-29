$machines{'vulcain'} = 
{ #soumission en interactif sur vulcain -> pas de 'submit'
    'nbcores'    => 8,
    'mempernode' => 32,
    #'execcmd'    => '',
    #'template'   => '',
    #'script'     => '',
    #'submit'     => '',
    #'time'       => "",
    #'args'       => "",
    #'argsbe'     => "",
    'bits'       => "_32bits",
    'hostarch'   => "i686_pc_linux",
    'ccprog'     => "gcc -Wall",
    'f77prog'    => "gfortran",
    'mpccprog'   => "mpicc -Wall",
    'mcfprog'    => "mpif90",
    'ccfopt'     => "-O3",
    'ccfdeb'     => "-O0 -g3",
    'f90flags'   => "-ffree-form -x f95-cpp-input",
    'ldflags'    => "-L/opt/intel/Compiler/11.0/083/mkl/lib/em64t/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lgfortran -lm -lrt",
    'ccflags'    => "",
    'lkfopt'     => "-s",
    'arprog'     => "ar",
    'arflags'    => "-ruv",
    'makecmd'    => "make"
}

