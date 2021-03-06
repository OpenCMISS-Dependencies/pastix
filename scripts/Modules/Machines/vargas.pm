$machines{'vargas'} =
{
    'nbcores'    => 32,
    'mempernode' => 102.4,
    'execcmd'    => '',
    'template'   => 'vargas.tpl',
    'script'     => 'job.par',
    'submit'     => 'llsubmit',
    'time'       => "01:00:00",
    'args'       => "",
    'argsbe'     => "",
    'bits'       => "_64bits",
    'hostarch'   => "power_ibm_aix",
    'ccprog'     => "xlc_r -ma -q64 -qlanglvl=extended -qarch=auto",
    'f77prog'    => "xlf_r -ma -q64 -qlanglvl=extended -qarch=auto",
    'mpccprog'   => "mpcc_r -ma -q64 -qlanglvl=extended -qarch=auto",
    'mcfprog'    => "mpxlf90_r",
    'ccfopt'     => " -O2 -qstrict -qtune=auto -s ",
    'ccfdeb'     => " -ma -q64 -qlanglvl=extended -qarch=auto ",
    'f90flags'   => "-qsuffix=cpp=f90 ",
    'ldflags'    => "-lessl -lxlf90 -lm -lrt ",
    'ccflags'    => "",
    'lkfopt'     => "-s",
    'arprog'     => "ar -X32_64 ",
    'arflags'    => "-ruv",
    'makecmd'    => "gmake"
}
