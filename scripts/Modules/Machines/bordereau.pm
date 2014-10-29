$machines{'bordereau'} = 
{
    'nbcores'    => 4,
    #'mempernode' => 32,
    #'execcmd'    => '',
    #'template'   => '',
    #'script'     => '',
    #'submit'     => '',
    #'time'       => "",
    #'args'       => "",
    #'argsbe'     => "",
    'bits'       => "_64bits",
    'hostarch'   => "i686_pc_linux",
    'ccprog'     => "gcc -Wall",
    'f77prog'    => "gfortran",
    'mpccprog'   => "mpicc -Wall",
    'mcfprog'    => "mpif90",
    'ccfopt'     => "-O3",
    'ccfdeb'     => "-O0 -g3",
    'f90flags'   => "-ffree-form -x f95-cpp-input",
    'ldflags'    => "-L/softs/acml/gnu64/lib -lacml -lm -lrt -lg2c",
    'ccflags'    => "",
    'lkfopt'     => "-s",
    'arprog'     => "ar",
    'arflags'    => "-ruv",
    'makecmd'    => "make"

}
