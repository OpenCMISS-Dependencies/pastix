    Program: PaStiX
    ___________________________________________________________________________

    About: Overview

	PaStiX (Parallel Sparse matriX package) is a scientific library that
    	provides a high performance parallel solver for very large sparse
    	linear systems based on direct methods. Numerical algorithms
	are implemented in single or double precision (real or complex)
	using LLt, LDLt and LU with static pivoting (for non symmetric
	matrices having a symmetric pattern). This solver provides also
	an adaptive blockwise iLU(k) factorization that can be used as a
    	parallel preconditioner using approximated supernodes to build
	a coarser block structure of the incomplete factors.

    About: Links

	HomePage of the project : http://pastix.gforge.inria.fr/

    	Downloads : https://gforge.inria.fr/frs/?group_id=186

   	Quick reference card :
	https://gforge.inria.fr/docman/index.php?group_id=186

    	Publications related to this project can be found on :
	http://www.labri.fr/perso/ramet/bib/Keyword/SPARSE.html

	Murge common solver interface :
	http://murge.gforge.inria.fr/


    Topic: Building library

	To build PaStiX library:
	     - copy one config file from src/config/ to config.in,
	     - edit this file to set your options,
	     - run
	       > make

	Check that your installation is working:
	     - run
	       > make examples
	       > ./example/bin/simple -lap 100

	After that, PaStiX can be used by:
	  - including :
	    > install/bin/pastix-conf --incs
	  - linking with :
	    > install/bin/pastix-conf --libs
	    or
	    > install/bin/pastix-conf --libs_murge

	See more information in the quick reference card.
	(https://gforge.inria.fr/docman/index.php?group_id=186)

    Topic: Calling PaStiX from C

	To call PaStiX from C refer to <pastix> function definition.

	You can get information about some PaStiX parameters
	by looking at <api.h>.

	You can also have a look at the examples in example/src directory.

    Topic: Calling PaStiX from Fortran

	To call PaStiX from Fortran refer to <pastix_fortran>.

	You can also have a look at the examples in example/src directory.

    Topic: PaStiX Team

	Please use this forum for submitting your problems :
	https://gforge.inria.fr/forum/?group_id=186

	To contact the PaStiX team :

	Casadei Astrid  : astrid.casadei@inria.fr

	Faverge Mathieu : faverge@labri.fr

	Lacoste Xavier  : xavier.lacoste@inria.fr

	Ramet Pierre    : ramet@labri.fr

    Topic: Some other solvers

	Direct solvers:

    	   - MUMPS : http://mumps.enseeiht.fr/
   	   - SuperLU : http://crd-legacy.lbl.gov/~xiaoye/SuperLU/
    	   - UMFPACK : http://www.cise.ufl.edu/research/sparse/umfpack/
    	   - TAUCS : http://www.tau.ac.il/~stoledo/taucs/
    	   - PARDISO : http://www.pardiso-project.org/
    	   - WSMP : http://researcher.ibm.com/view_project.php?id=1426
    	   - PSPASES : http://www-users.cs.umn.edu/~mjoshi/pspases/
    	   - HSL : http://www.hsl.rl.ac.uk/
    	   - BCSLib : http://www.boeing.com/phantom/bcslib/

	Hybrid solvers:

    	   - HIPS : http://hips.gforge.inria.fr/
    	   - pARMS : http://www-users.cs.umn.edu/~saad/software/pARMS/
    	   - MaPHYS : (no website yet)
    	   - PDSLin : (no website yet)

	Generic toolkits:

    	   - PETSc : http://www.mcs.anl.gov/petsc/
    	   - Trilinos : http://trilinos.sandia.gov/

    About: License

       PaStiX is distributed under the CeCILL-C license.

       > Copyright 2008 BORDEAUX I UNIVERSITY & INRIA
       >
       > This software is governed by the CeCILL-C license under French law
       > and abiding by the rules of distribution of free software. You can
       > use, modify and/or redistribute the software under the terms of the
       > CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
       > URL: "http://www.cecill.info".
       >
       > As a counterpart to the access to the source code and rights to copy,
       > modify and redistribute granted by the license, users are provided
       > only with a limited warranty and the software's author, the holder of
       > the economic rights, and the successive licensors have only limited
       > liability.
       >
       > In this respect, the user's attention is drawn to the risks associated
       > with loading, using, modifying and/or developing or reproducing the
       > software by the user in light of its specific status of free software,
       > that may mean that it is complicated to manipulate, and that also
       > therefore means that it is reserved for developers and experienced
       > professionals having in-depth computer knowledge. Users are therefore
       > encouraged to load and test the software's suitability as regards
       > their requirements in conditions enabling the security of their
       > systems and/or data to be ensured and, more generally, to use and
       > operate it in the same conditions as regards security.
       >
       > The fact that you are presently reading this means that you have had
       > knowledge of the CeCILL-C license and that you accept its terms.
