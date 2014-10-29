$machines{'argent'} = 
{
    'nbcores'    => 8,
    'mempernode' => 0,
    'execcmd'    => '',
    'template'   => 'argent.tpl',
    'script'     => 'job.par',
    'submit'     => 'bsub <',
    'time'       => 60,
    'args'       => "",
    'argsbe'     => "",
    'bits'       => "_64bits",
    'hostarch'   => "i686_pc_linux",
    'ccprog'     => "icc -Wall",
    'f77prog'    => "ifort",
    'mpccprog'   => "mpicc -Wall",
    'mcfprog'    => "mpif90",
    'ccfopt'     => "-O3",
    'ccfdeb'     => "-O0 -g3",
    'f90flags'   => "-fpp",
    'ldflags'    => "-L/applications/intel/ict/3.0/cmkl/9.0/lib/64 -lmkl -L/opt/intelruntime -lifcore -lm -lrt ",
    'ccflags'    => " -DPASTIX_FUNNELED",
    'lkfopt'     => "",
    'arprog'     => "ar ",
    'arflags'    => "-ruv",
    'makecmd'    => "make"

}
