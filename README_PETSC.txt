Topic: To use PaStiX through PETSc

       Configure PETSc :
       > export PETSC_DIR=$PWD
       > ./configure --download-pastix=1 --download-ptscotch=1

       Building PETSc :
       > export PETSC_DIR=$PWD
       > make all test

       Testing PaStiX through PETSc :
       > cd src/ksp/pc/examples/tests/
       > make ex2
       > ./ex2 -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package pastix -ksp_monitor \
       >       -mat_pastix_threadnbr [1:X]  -mat_pastix_verbose [0:4] -mat_pastix_check [0,1]
