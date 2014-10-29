c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     routines from SPARSKIT and extensions 
c-----------------------------------------------------------------------
c----------------------------------------------------------------------- 
      subroutine dperm (nrow,a,ja,ia,ao,jao,iao,perm,qperm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),
     +        qperm(*),job
      real*8 a(*),ao(*) 
c-----------------------------------------------------------------------
c This routine permutes the rows and columns of a matrix stored in CSR
c format. i.e., it computes P A Q, where P, Q are permutation matrices. 
c P maps row i into row perm(i) and Q maps column j into column qperm(j): 
c      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix
c In the particular case where Q is the transpose of P (symmetric 
c permutation of A) then qperm is not needed. 
c note that qperm should be of length ncol (number of columns) but this
c is not checked. 
c-----------------------------------------------------------------------
c Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991. 
c-----------------------------------------------------------------------
c on entry: 
c---------- 
c n 	= dimension of the matrix
c a, ja, 
c    ia = input matrix in a, ja, ia format
c perm 	= integer array of length n containing the permutation arrays
c	  for the rows: perm(i) is the destination of row i in the
c         permuted matrix -- also the destination of column i in case
c         permutation is symmetric (job .le. 2) 
c
c qperm	= same thing for the columns. This should be provided only
c         if job=3 or job=4, i.e., only in the case of a nonsymmetric
c	  permutation of rows and columns. Otherwise qperm is a dummy
c
c job	= integer indicating the work to be done:
c * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P)
c 		job = 1	permute a, ja, ia into ao, jao, iao 
c 		job = 2 permute matrix ignoring real values.
c * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q 
c 		job = 3	permute a, ja, ia into ao, jao, iao 
c 		job = 4 permute matrix ignoring real values.
c		
c on return: 
c-----------
c ao, jao, iao = input matrix in a, ja, ia format
c
c in case job .eq. 2 or job .eq. 4, a and ao are never referred to 
c and can be dummy arguments. 
c Notes:
c------- 
c  1) algorithm is in place 
c  2) column indices may not be sorted on return even  though they may be 
c     on entry.
c----------------------------------------------------------------------c
c local variables 
      integer locjob, mod
c
c     locjob indicates whether or not real values must be copied. 
c     
      locjob = mod(job,2) 
c
c permute rows first 
c 
      call rperm (nrow,a,ja,ia,ao,jao,iao,perm,locjob)
c
c then permute columns
c
      locjob = 0
c
      if (job .le. 2) then
         call cperm (nrow,ao,jao,iao,ao,jao,iao,perm,locjob) 
      else 
         call cperm (nrow,ao,jao,iao,ao,jao,iao,qperm,locjob) 
      endif 
c     
      return
c-------end-of-dperm----------------------------------------------------
c-----------------------------------------------------------------------
      end

c-----------------------------------------------------------------------
      subroutine dperm1 (i1,i2,a,ja,ia,b,jb,ib,perm,ipos,job)
      integer i1,i2,job,ja(*),ia(*),jb(*),ib(*),perm(*)
      real*8 a(*),b(*)
c----------------------------------------------------------------------- 
c     general submatrix permutation/ extraction routine.
c----------------------------------------------------------------------- 
c     extracts rows perm(i1), perm(i1+1), ..., perm(i2) (in this order) 
c     from a matrix (doing nothing in the column indices.) The resulting 
c     submatrix is constructed in b, jb, ib. A pointer ipos to the
c     beginning of arrays b,jb,is also allowed (i.e., nonzero elements
c     are accumulated starting in position ipos of b, jb). 
c-----------------------------------------------------------------------
c Y. Saad,Sep. 21 1989 / recoded Jan. 28 1991 / modified for SPARSLIB 
c Sept. 1997.. 
c-----------------------------------------------------------------------
c on entry: 
c---------- 
c n 	= dimension of the matrix
c a,ja,
c   ia  = input matrix in CSR format
c perm 	= integer array of length n containing the indices of the rows
c         to be extracted. 
c
c job   = job indicator. if (job .ne.1) values are not copied (i.e.,
c         only pattern is copied).
c
c on return: 
c-----------
c b,jb,
c ib   = matrix in csr format. b(ipos:ipos+nnz-1),jb(ipos:ipos+nnz-1) 
c     contain the value and column indices respectively of the nnz
c     nonzero elements of the permuted matrix. thus ib(1)=ipos.
c
c Notes:
c------- 
c  algorithm is NOT in place 
c----------------------------------------------------------------------- 
c local variables
c
      integer ko,irow,k 
      logical values
c-----------------------------------------------------------------------
      values = (job .eq. 1) 
      ko = ipos 
      ib(1) = ko
      do 900 i=i1,i2
         irow = perm(i) 
         do 800 k=ia(irow),ia(irow+1)-1
            if (values) b(ko) = a(k)
            jb(ko) = ja(k)
            ko=ko+1
 800     continue
         ib(i-i1+2) = ko
 900  continue
      return
c--------end-of-dperm1--------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine dperm2 (i1,i2,a,ja,ia,b,jb,ib,perm,rperm,istart,
     *        ipos,job)
      integer i1,i2,job,istart,ja(*),ia(*),jb(*),ib(*),perm(*),rperm(*) 
      real*8 a(*),b(*)
c----------------------------------------------------------------------- 
c     general submatrix permutation/ extraction routine.
c----------------------------------------------------------------------- 
c     extracts rows rperm(i1), rperm(i1+1), ..., rperm(i2) and does an
c     associated column permutation (using array perm). The resulting
c     submatrix is constructed in b, jb, ib. For added flexibility,
c     the extracted are put in sequence starting from row 'istart' of B. 
c     In addition a pointer ipos to the beginning of arrays b,jb,is also 
c     allowed (i.e., nonzero elements are accumulated starting in
c     position ipos of b, jb). 
c     extremely flexible. exple
c     to permute msr to msr (excluding diagonals) 
c     call dperm2 (1,n,a,ja,ja,b,jb,jb,perm,rperm,1,n+2) 
c-----------------------------------------------------------------------
c Y. Saad,Sep. 21 1989 / recoded Jan. 28 1991. 
c-----------------------------------------------------------------------
c on entry: 
c---------- 
c n 	= dimension of the matrix
c a,ja,
c   ia  = input matrix in CSR format
c perm 	= integer array of length n containing the permutation arrays
c	  for the rows: perm(i) is the destination of row i in the
c         permuted matrix -- also the destination of column i. 
c
c rperm	= reverse permutation. defined by rperm(iperm(i))=i for all i
c
c job   = job indicator. if (job .ne.1) values are not copied (i.e.,
c         only pattern is copied).
c
c on return: 
c-----------
c b,ja,
c ib   = matrix in csr format. positions 1,2,...,istart-1 of ib 
c     are not touched. b(ipos:ipos+nnz-1),jb(ipos:ipos+nnz-1) 
c     contain the value and column indices respectively of the nnz
c     nonzero elements of the permuted matrix. thus ib(istart)=ipos.
c
c Notes:
c------- 
c  1) algorithm is NOT in place 
c  2) column indices may not be sorted on return even  though they 
c     may be on entry.
c----------------------------------------------------------------------- 
c local variables
c
      integer ko,irow,k 
      logical values
c-----------------------------------------------------------------------
      values = (job .eq. 1) 
      ko = ipos 
      ib(istart) = ko
      do 900 i=i1,i2
         irow = rperm(i) 
         do 800 k=ia(irow),ia(irow+1)-1
            if (values) b(ko) = a(k)
            jb(ko) = perm(ja(k))
            ko=ko+1
 800     continue
         ib(istart+i-i1+1) = ko
 900  continue
      return
c--------end-of-dperm2--------------------------------------------------
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine multic (n,ja,ia,ncol,kolrs,il,iord,maxcol,ierr) 
      integer n, ja(*),ia(n+1),kolrs(n),iord(n),il(maxcol+1),ierr 
c-----------------------------------------------------------------------
c multic:
c     multicoloring ordering -- greedy algorithm -- 
c     determines the coloring permutation and sets up
c     corresponding data structures for it.
c-----------------------------------------------------------------------
c on entry
c -------- 
c n     = row and column dimention of matrix
c ja    = column indices of nonzero elements of matrix, stored rowwise.
c ia    = pointer to beginning of each row in ja.
c maxcol= maximum number of colors allowed -- the size of il is
c         maxcol+1 at least. Note: the number of colors does not
c         exceed the maximum degree of each node +1.
c 
c on return
c --------- 
c ncol  = number of colours found 
c kolrs = integer array containing the color number assigned to each node 
c il    = integer array containing the pointers to the
c         beginning of each color set. In the permuted matrix
c         the rows /columns il(kol) to il(kol+1)-1 have the same color.
c iord  = permutation array corresponding to the multicolor ordering.
c         row number i will become row nbumber iord(i) in permuted 
c         matrix. (iord = destination permutation array).
c ierr  = integer. Error message. normal return ierr = 0
c-----------------------------------------------------------------------
c     
      integer kol, i, j, k, maxcol, mycol 
c     
      ierr = 0
      do 1 j=1, n
         kolrs(j) = 0
         iord(j) = j
 1    continue 
      do 11 j=1, maxcol
         il(j) = 0
 11   continue
c     
      ncol = 0
c     
c     scan all nodes 
c     
      do 4 ii=1, n
         i = iord(ii) 
c     
c     look at adjacent nodes to determine colors already assigned
c     
         mcol = 0
         do 2 k=ia(i), ia(i+1)-1
            j = ja(k)
            icol = kolrs(j)
            if (icol .ne. 0) then
               mcol = max(mcol,icol) 
c     
c     il used as temporary to record already assigned colors.
c     
               il(icol) = 1 
            endif
 2       continue
c     
c     colors taken determined. scan il until a slot opens up.
c     
         mycol = 1
 3       if (il(mycol) .eq. 1) then
            mycol = mycol+1 
c
c #### put test elsewhere --
c 
            if (mycol .gt. maxcol) goto 99
            if (mycol .le. mcol) goto 3
         endif
c     
c     reset il to zero for next nodes
c     
         do 35 j=1, mcol
            il(j) = 0
 35      continue
c     
c     assign color and update number of colors so far
c     
         kolrs(i) = mycol
         ncol = max(ncol,mycol)
 4    continue
c     
c     every node has now been colored. Count nodes of each color
c     
      do 6 j=1, n
         kol = kolrs(j)+1
         il(kol) = il(kol)+1
 6    continue
c     
c     set pointers il
c     
      il(1) = 1
      do 7 j=1, ncol
         il(j+1) = il(j)+il(j+1)
 7    continue
c     
c     set iord
c     
      do 8 j=1, n
         kol = kolrs(j) 
         iord(j) = il(kol)
         il(kol) = il(kol)+1
 8    continue
c     
c     shift il back 
c     
      do 9 j=ncol,1,-1
         il(j+1) = il(j)
 9    continue
      il(1) = 1
c     
      return
 99   ierr = 1
      return
c-----end-of-multic-----------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine xtrows (i1,i2,a,ja,ia,ao,jao,iao,iperm,job)
      integer i1,i2,ja(*),ia(*),jao(*),iao(*),iperm(*),job
      real*8 a(*),ao(*) 
c-----------------------------------------------------------------------
c this subroutine extracts given rows from a matrix in CSR format. 
c Specifically, rows number iperm(i1), iperm(i1+1), ...., iperm(i2)
c are extracted and put in the output matrix ao, jao, iao, in CSR
c format.  NOT in place. 
c Youcef Saad -- coded Feb 15, 1992. 
c-----------------------------------------------------------------------
c on entry:
c----------
c i1,i2   = two integers indicating the rows to be extracted.
c           xtrows will extract rows iperm(i1), iperm(i1+1),..,iperm(i2),
c           from original matrix and stack them in output matrix 
c           ao, jao, iao in csr format
c
c a, ja, ia = input matrix in csr format
c
c iperm	= integer array of length nrow containing the reverse permutation 
c         array for the rows. row number iperm(j) in permuted matrix PA
c         used to be row number j in unpermuted matrix.
c         ---> a(i,j) in the permuted matrix was a(iperm(i),j) 
c         in the inout matrix.
c
c job	= integer indicating the work to be done:
c 		job .ne. 1 : get structure only of output matrix,,
c               i.e., ignore real values. (in which case arrays a 
c               and ao are not used nor accessed).
c 		job = 1	get complete data structure of output matrix. 
c               (i.e., including arrays ao and iao).
c------------
c on return: 
c------------ 
c ao, jao, iao = input matrix in a, ja, ia format
c note : 
c        if (job.ne.1)  then the arrays a and ao are not used.
c----------------------------------------------------------------------c
c           Y. Saad, revised May  2, 1990                              c
c----------------------------------------------------------------------c
      logical values
      values = (job .eq. 1) 
c
c copying 
c
      ko = 1
      iao(1) = ko
      do 100 j=i1,i2 
c
c ii=iperm(j) is the index of old row to be copied.
c        
         ii = iperm(j) 
         do 60 k=ia(ii), ia(ii+1)-1 
            jao(ko) = ja(k) 
            if (values) ao(ko) = a(k)
            ko = ko+1
 60      continue
         iao(j-i1+2) = ko
 100  continue
c
      return
c---------end-of-xtrows------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      subroutine amub (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *                 c,jc,ic,nzmax,iw,ierr) 
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(ncol+1),
     *     ic(ncol+1),iw(ncol)
c-----------------------------------------------------------------------
c     performs the matrix by matrix product C = A B 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow  = integer. The row dimension of A
c ncol  = integer. The column dimension of A
c job   = integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c b, 
c jb, 
c ib    =  Matrix B in compressed sparse row format.
c
c nzmax = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic    = resulting matrix C in compressed sparse row sparse format.
c           
c ierr  = integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw    = integer work array of length equal to the number of
c         columns in A.
c Notes:
c-------
c         The column dimension of B is not needed.
c
c-----------------------------------------------------------------
      real*8 scal 
      logical values
      values = (job .ne. 0) 
      len = 0
      ic(1) = 1 
      ierr = 0
c     initialize array iw.
      do 1 j=1, ncol
         iw(j) = 0
 1           continue
c
      do 500 ii=1, nrow 
c     row i 
         do 200 ka=ia(ii), ia(ii+1)-1 
            if (values) scal = a(ka)
            jj   = ja(ka)
            do 100 kb=ib(jj),ib(jj+1)-1
               jcol = jb(kb)
               jpos = iw(jcol)
               if (jpos .eq. 0) then
                  len = len+1
                  if (len .gt. nzmax) then
                     ierr = ii
                     return
                  endif
                  jc(len) = jcol
                  iw(jcol)= len
                  if (values) c(len)  = scal*b(kb)
               else
                  if (values) c(jpos) = c(jpos) + scal*b(kb)
               endif
 100                   continue
 200                        continue
         do 201 k=ic(ii), len
            iw(jc(k)) = 0
 201             continue
         ic(ii+1) = len+1
 500       continue
      return
c-------------end-of-amub-----------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine dump (i1,i2,values,a,ja,ia,iout)
      integer i1, i2, ia(*), ja(*), iout
      real*8 a(*) 
      logical values
c outputs rows i1 through i2 of a sparse matrix stored in CSR format 
c (or columns i1 through i2 of a matrix stored in CSC format) in a file, 
c one (column) row at a time in a nice readable format. 
c This is a simple routine which is useful for debugging. 
c-----------------------------------------------------------------------
c on entry:
c---------
c i1    = first row (column) to print out
c i2    = last row (column) to print out 
c values= logical. indicates whether or not to print real values.
c         if value = .false. only the pattern will be output.
c a,
c ja, 
c ia    =  matrix in CSR format (or CSC format) 
c iout  = logical unit number for output.
c---------- 
c the output file iout will have written in it the rows or columns 
c of the matrix in one of two possible formats (depending on the max 
c number of elements per row. The values are output with only 
c two digits of accuracy (D9.2). )
c-----------------------------------------------------------------------
c     local variables
      integer maxr, i, k1, k2 
c
c select mode horizontal or vertical 
c
        maxr = 0
        do 1 i=i1, i2
           maxr = max0(maxr,ia(i+1)-ia(i))
 1               continue
        
        if (maxr .le. 8) then
c
c able to do one row acros line
c
        do 2 i=i1, i2
           write(iout,100) i
           k1=ia(i)
           k2 = ia(i+1)-1
           write (iout,101) (ja(k),k=k1,k2)
           if (values) write (iout,102) (a(k),k=k1,k2)
 2               continue
      else 
c
c unable to one row acros line. do three items at a time
c across a line 
         do 3 i=i1, i2
            if (values) then
               write(iout,200) i
            else
               write(iout,203) i               
            endif
            k1=ia(i)
            k2 = ia(i+1)-1
            if (values) then
               write (iout,201) (ja(k),a(k),k=k1,k2)
            else
                write (iout,202) (ja(k),k=k1,k2)
             endif
 3                  continue
      endif 
c
c formats :
c
 100    format (1h ,35(1h-),' row',i6,1x,35(1h-) )
 101      format(' col:',8(i5,6h     : ))
 102        format(' val:',8(D9.2,2h :) )
 200          format (1h ,31(1h-),' row',i3,1x,31(1h-),/
     *     3('  columns :    values  * ') )
c-------------xiiiiiihhhhhhddddddddd-*-
 201            format(3(1h ,i6,6h   :  ,D9.2,3h * ) )
 202              format(6(1h ,i5,6h  *    ) ) 
 203                format (1h ,31(1h-),' row',i3,1x,31(1h-),/
     *     3('  column  :  column   *') )
      return
c----end-of-dump--------------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine readmt (nmax,nzmax,job,iounit,a,ja,ia,rhs,nrhs,
     *     guesol,nrow,ncol,nnz,title,key,type,ierr)
c-----------------------------------------------------------------------
c     this subroutine reads  a boeing/harwell matrix. handles right hand 
c     sides in full format only (no sparse right hand sides).
c     Also the matrix must be in assembled forms.
c     Author: Youcef Saad - Date: Sept., 1989
c     updated Oct 31, 1989.
c-----------------------------------------------------------------------
c on entry:
c---------
c nmax 	 =  max column dimension  allowed for matrix. The array ia should 
c	    be of length at least ncol+1 (see below) if job.gt.0
c nzmax	 = max number of nonzeros elements allowed. the arrays a, 
c          and ja should be of length equal to nnz (see below) if these
c          arrays are to be read (see job).
c          
c job	 = integer to indicate what is to be read. (note: job is an
c          input and output parameter, it can be modified on return)
c          job = 0    read the values of ncol, nrow, nnz, title, key,
c                     type and return. matrix is not read and arrays
c                     a, ja, ia, rhs are not touched.
c          job = 1    read srtucture only, i.e., the arrays ja and ia.
c          job = 2    read matrix including values, i.e., a, ja, ia
c          job = 3    read matrix and right hand sides: a,ja,ia,rhs.
c		      rhs may contain initial guesses and exact 
c                     solutions appended to the actual right hand sides.
c		      this will be indicated by the output parameter
c                     guesol [see below]. 
c                     
c nrhs   = integer. nrhs is an input as well as ouput parameter.
c          at input nrhs contains the total length of the array rhs.
c          See also ierr and nrhs in output parameters.
c
c iounit = logical unit number where to read the matrix from.
c
c on return:
c---------- 
c job    = on return job may be modified to the highest job it could
c          do: if job=2 on entry but no matrix values are available it
c          is reset to job=1 on return. Similarly of job=3 but no rhs 
c          is provided then it is rest to job=2 or job=1 depending on 
c          whether or not matrix values are provided.
c          Note that no error message is triggered (i.e. ierr = 0 
c          on return in these cases. It is therefore important to
c          compare the values of job on entry and return ).
c
c a	 = the a matrix in the a, ia, ja (column) storage format
c ja 	 = row number of element a(i,j) in array a.
c ia     = pointer  array. ia(i) points to the beginning of column i.
c
c rhs    = real array of size nrow + 1 if available (see job)
c
c nrhs   = integer containing the number of right-hand sides found
c          each right hand side may be accompanied with an intial guess
c          and also the exact solution.
c
c guesol = a 2-character string indicating whether an initial guess 
c          (1-st character) and / or the exact solution (2-nd
c          character) is provided with the right hand side.
c	   if the first character of guesol is 'G' it means that an
c          an intial guess is provided for each right-hand side.
c          These are appended to the right hand-sides in the array rhs.
c	   if the second character of guesol is 'X' it means that an
c          exact solution is provided for each right-hand side.
c          These are  appended to the right hand-sides 
c          and the initial guesses (if any) in the array rhs.
c
c nrow   = number of rows in matrix
c ncol	 = number of columns in matrix 
c nnz	 = number of nonzero elements in A. This info is returned
c          even if there is not enough space in a, ja, ia, in order
c          to determine the minimum storage needed. 
c
c title  = character*100 = title of matrix test ( character a*72). 
c key    = character*8  = key of matrix 
c type   = charatcer*3  = type of matrix.
c          for meaning of title, key and type refer to documentation 
c          Harwell/Boeing matrices.
c
c ierr   = integer used for error messages 
c         * ierr  =  0 means that  the matrix has been read normally. 
c         * ierr  =  1 means that  the array matrix could not be read 
c         because ncol+1 .gt. nmax
c         * ierr  =  2 means that  the array matrix could not be read 
c         because nnz .gt. nzmax 
c         * ierr  =  3 means that  the array matrix could not be read 
c         because both (ncol+1 .gt. nmax) and  (nnz .gt. nzmax )
c         * ierr  =  4 means that  the right hand side (s) initial 
c         guesse (s) and exact solution (s)   could  not be
c         read because they are stored in sparse format (not handled
c         by this routine ...) 
c         * ierr  =  5 means that the right-hand-sides, initial guesses
c         and exact solutions could not be read because the length of 
c         rhs as specified by the input value of nrhs is not 
c         sufficient to store them. The rest of the matrix may have
c         been read normally.
c 
c Notes:
c-------
c 1) The file inout must be open (and possibly rewound if necessary)
c    prior to calling readmt.
c 2) Refer to the documentation on the Harwell-Boeing formats
c    for details on the format assumed by readmt.
c    We summarize the format here for convenience.
c  
c    a) all lines in inout are assumed to be 80 character long.
c    b) the file consists of a header followed by the block of the 
c       column start pointers followed by the block of the
c       row indices, followed by the block of the real values and
c       finally the numerical values of the right-hand-side if a 
c       right hand side is supplied. 
c    c) the file starts by a header which contains four lines if no
c       right hand side is supplied and five lines otherwise.
c       * first line contains the title (72 characters long) followed by
c         the 8-character identifier (name of the matrix, called key)
c        [ A72,A8 ]
c       * second line contains the number of lines for each
c         of the following data blocks (4 of them) and the total number 
c         of lines excluding the header.
c        [5i4]
c       * the third line contains a three character string identifying
c         the type of matrices as they are referenced in the Harwell
c         Boeing documentation [e.g., rua, rsa,..] and the number of
c         rows, columns, nonzero entries.
c         [A3,11X,4I14]
c       * The fourth line contains the variable fortran format
c         for the following data blocks.
c         [2A16,2A20] 
c       * The fifth line is only present if right-hand-sides are 
c         supplied. It consists of three one character-strings containing
c         the storage format for the right-hand-sides 
c         ('F'= full,'M'=sparse=same as matrix), an initial guess 
c         indicator ('G' for yes), an exact solution indicator 
c         ('X' for yes), followed by the number of right-hand-sides
c         and then the number of row indices. 
c         [A3,11X,2I14] 
c     d) The three following blocks follow the header as described 
c        above.
c     e) In case the right hand-side are in sparse formats then 
c        the fourth block uses the same storage format as for the matrix
c        to describe the NRHS right hand sides provided, with a column
c        being replaced by a right hand side.
c-----------------------------------------------------------------------
        character title*72, key*8, type*3, ptrfmt*16, indfmt*16,
     1     valfmt*20, rhsfmt*20, rhstyp*3, guesol*2
      integer totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol,
     1     nnz, neltvl, nrhs, nmax, nzmax, nrwindx
      integer ia (nmax+1), ja (nzmax) 
      real*8 a(nzmax), rhs(*) 
c-----------------------------------------------------------------------
      ierr = 0
      lenrhs = nrhs
c     
      read (iounit,10) title, key, totcrd, ptrcrd, indcrd, valcrd, 
     1     rhscrd, type, nrow, ncol, nnz, neltvl, ptrfmt, indfmt, 
     2     valfmt, rhsfmt
 10   format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
c     
      if (rhscrd .gt. 0) read (iounit,11) rhstyp, nrhs, nrwindx
 11   format (a3,11x,i14,i14)
c     
c     anything else to read ?
c     
      if (job .le. 0) return
c     ---- check whether matrix is readable ------ 
      n = ncol
      if (ncol .gt. nmax) ierr = 1
      if (nnz .gt. nzmax) ierr = ierr + 2
      if (ierr .ne. 0) return
c     ---- read pointer and row numbers ---------- 
      read (iounit,ptrfmt) (ia (i), i = 1, n+1)
      read (iounit,indfmt) (ja (i), i = 1, nnz)
c     --- reading values of matrix if required....
      if (job .le. 1)  return
c     --- and if available ----------------------- 
      if (valcrd .le. 0) then
	 job = 1
	 return
      endif
      read (iounit,valfmt) (a(i), i = 1, nnz)
c     --- reading rhs if required ---------------- 
      if (job .le. 2)  return
c     --- and if available ----------------------- 
      if ( rhscrd .le. 0) then
	 job = 2
	 return
      endif
c     
c     --- read right-hand-side.-------------------- 
c     
      if (rhstyp(1:1) .eq. 'M') then 
         ierr = 4
         return
      endif
c     
      guesol = rhstyp(2:3) 
c     
      nvec = 1 
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') nvec=nvec+1
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') nvec=nvec+1
c     
      len = nrhs*nrow 
c     
      if (len*nvec .gt. lenrhs) then
         ierr = 5
         return
      endif
c     
c     read right-hand-sides
c     
      next = 1
      iend = len
      read(iounit,rhsfmt) (rhs(i), i = next, iend)
c     
c     read initial guesses if available
c     
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') then
         next = next+len
         iend = iend+ len
         read(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c     
c     read exact solutions if available
c     
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') then
         next = next+len
         iend = iend+ len
         read(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c     
      return
c---------end-of-readmt-------------------------------------------------
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine roscal(nrow,job,nrm,a,ja,ia,diag,b,jb,ib,ierr) 
      real*8 a(*), b(*), diag(nrow) 
      integer nrow,job,nrm,ja(*),jb(*),ia(nrow+1),ib(nrow+1),ierr 
c-----------------------------------------------------------------------
c scales the rows of A such that their norms are one on return
c 3 choices of norms: 1-norm, 2-norm, max-norm.
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = diagonal matrix stored as a vector containing the matrix
c        by which the rows have been scaled, i.e., on return 
c        we have B = Diag*A.
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c	    
c ierr  = error message. ierr=0     : Normal return 
c                        ierr=i > 0 : Row number i is a zero row.
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)        algorithm in place (B can take the place of A).
c-----------------------------------------------------------------
      call rnrms (nrow,nrm,a,ja,ia,diag)

      ierr = 0
      do 1 j=1, nrow
         if (diag(j) .eq. 0.0d0) then
            ierr = j 
            return
         else
            diag(j) = 1.0d0/diag(j)
         endif
 1    continue
      call diamua(nrow,job,a,ja,ia,diag,b,jb,ib)
      return
c-------end-of-roscal---------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine coscal(nrow,job,nrm,a,ja,ia,diag,b,jb,ib,ierr) 
c----------------------------------------------------------------------- 
      real*8 a(*),b(*),diag(nrow) 
      integer nrow,job,ja(*),jb(*),ia(nrow+1),ib(nrow+1),ierr 
c-----------------------------------------------------------------------
c scales the columns of A such that their norms are one on return
c result matrix written on b, or overwritten on A.
c 3 choices of norms: 1-norm, 2-norm, max-norm. in place.
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = diagonal matrix stored as a vector containing the matrix
c        by which the columns have been scaled, i.e., on return 
c        we have B = A * Diag
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c
c ierr  = error message. ierr=0     : Normal return 
c                        ierr=i > 0 : Column number i is a zero row.
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)       algorithm in place (B can take the place of A).
c-----------------------------------------------------------------
      call cnrms (nrow,nrm,a,ja,ia,diag)
      ierr = 0
      do 1 j=1, nrow
         if (diag(j) .eq. 0.0) then
            ierr = j 
            return
         else
            diag(j) = 1.0d0/diag(j)
         endif
 1    continue
      call amudia (nrow,job,a,ja,ia,diag,b,jb,ib)
      return
c--------end-of-coscal-------------------------------------------------- 
c-----------------------------------------------------------------------
      end
      subroutine rnrms   (nrow, nrm, a, ja, ia, diag) 
      real*8 a(*), diag(nrow), scal 
      integer ja(*), ia(nrow+1) 
c-----------------------------------------------------------------------
c gets the norms of each row of A. (choice of three norms)
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = real vector of length nrow containing the norms
c
c-----------------------------------------------------------------
      do 1 ii=1,nrow
c
c     compute the norm if each element.
c     
         scal = 0.0d0
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         if (nrm .eq. 0) then
            do 2 k=k1, k2
               scal = max(scal,abs(a(k) ) ) 
 2          continue
         elseif (nrm .eq. 1) then
            do 3 k=k1, k2
               scal = scal + abs(a(k) ) 
 3          continue
         else
            do 4 k=k1, k2
               scal = scal+a(k)**2
 4          continue
         endif 
         if (nrm .eq. 2) scal = sqrt(scal) 
         diag(ii) = scal
 1    continue
      return
c-----------------------------------------------------------------------
c-------------end-of-rnrms----------------------------------------------
      end 
c-----------------------------------------------------------------------
      subroutine cnrms   (nrow, nrm, a, ja, ia, diag) 
      real*8 a(*), diag(nrow) 
      integer ja(*), ia(nrow+1) 
c-----------------------------------------------------------------------
c gets the norms of each column of A. (choice of three norms)
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c on return:
c----------
c
c diag = real vector of length nrow containing the norms
c
c-----------------------------------------------------------------
      do 10 k=1, nrow 
         diag(k) = 0.0d0
 10   continue
      do 1 ii=1,nrow
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         do 2 k=k1, k2
            j = ja(k) 
c     update the norm of each column
            if (nrm .eq. 0) then
               diag(j) = max(diag(j),abs(a(k) ) ) 
            elseif (nrm .eq. 1) then
               diag(j) = diag(j) + abs(a(k) ) 
            else
               diag(j) = diag(j)+a(k)**2
            endif 
 2       continue
 1    continue
      if (nrm .ne. 2) return
      do 3 k=1, nrow
         diag(k) = sqrt(diag(k))
 3    continue
      return
c-----------------------------------------------------------------------
c------------end-of-cnrms-----------------------------------------------
      end 
      subroutine diamua (nrow,job, a, ja, ia, diag, b, jb, ib)
      real*8 a(*), b(*), diag(nrow), scal
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1) 
c-----------------------------------------------------------------------
c performs the matrix by matrix product B = Diag * A  (in place) 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c----------
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c	    
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)        algorithm in place (B can take the place of A).
c           in this case use job=0.
c-----------------------------------------------------------------
      do 1 ii=1,nrow
c     
c     normalize each row 
c     
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         scal = diag(ii) 
         do 2 k=k1, k2
            b(k) = a(k)*scal
 2       continue
 1    continue
c     
      if (job .eq. 0) return
c     
      do 3 ii=1, nrow+1
         ib(ii) = ia(ii)
 3    continue
      do 31 k=ia(1), ia(nrow+1) -1 
         jb(k) = ja(k)
 31   continue
      return
c----------end-of-diamua------------------------------------------------
c-----------------------------------------------------------------------
      end 
c----------------------------------------------------------------------- 
      subroutine amudia (nrow,job, a, ja, ia, diag, b, jb, ib)
      real*8 a(*), b(*), diag(nrow) 
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1) 
c-----------------------------------------------------------------------
c performs the matrix by matrix product B = A * Diag  (in place) 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c----------
c
c b, 
c jb, 
c ib	= resulting matrix B in compressed sparse row sparse format.
c	    
c Notes:
c-------
c 1)        The column dimension of A is not needed. 
c 2)        algorithm in place (B can take the place of A).
c-----------------------------------------------------------------
      do 1 ii=1,nrow
c     
c     scale each element 
c     
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         do 2 k=k1, k2
            b(k) = a(k)*diag(ja(k)) 
 2       continue
 1    continue
c     
      if (job .eq. 0) return
c     
      do 3 ii=1, nrow+1
         ib(ii) = ia(ii)
 3    continue
      do 31 k=ia(1), ia(nrow+1) -1 
         jb(k) = ja(k)
 31   continue
      return
c-----------------------------------------------------------------------
c-----------end-of-amudiag----------------------------------------------
      end 
c aplb   :   computes     C = A+B                                      c
      subroutine aplb (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1),
     *     iw(ncol)
c-----------------------------------------------------------------------
c performs the matrix sum  C = A+B. 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c job   = integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c b, 
c jb, 
c ib	=  Matrix B in compressed sparse row format.
c
c nzmax	= integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic	= resulting matrix C in compressed sparse row sparse format.
c	    
c ierr	= integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw	= integer work array of length equal to the number of
c         columns in A.
c
c-----------------------------------------------------------------------
      logical values

      values = (job .ne. 0) 
      ierr = 0
      len = 0
      ic(1) = 1 
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
c     
      do 500 ii=1, nrow
c     row i 
         do 200 ka=ia(ii), ia(ii+1)-1 
            len = len+1
            jcol    = ja(ka)
            if (len .gt. nzmax) goto 999
            jc(len) = jcol 
            if (values) c(len)  = a(ka) 
            iw(jcol)= len
 200     continue
c     
         do 300 kb=ib(ii),ib(ii+1)-1
            jcol = jb(kb)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
               len = len+1
               if (len .gt. nzmax) goto 999
               jc(len) = jcol
               if (values) c(len)  = b(kb)
               iw(jcol)= len
            else
               if (values) c(jpos) = c(jpos) + b(kb)
            endif
 300     continue
         do 301 k=ic(ii), len
	    iw(jc(k)) = 0
 301     continue
         ic(ii+1) = len+1
 500  continue
      return
 999  ierr = ii
      return
c------------end of aplb ----------------------------------------------- 
      end

c csrcsc  : converts compressed sparse row format to compressed sparse c
      subroutine csrcsc (n,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n+1),ja(*),jao(*)
      real*8  a(*),ao(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to      Compressed Sparse Column
c
c (transposition operation)   Not in place. 
c----------------------------------------------------------------------- 
c -- not in place --
c this subroutine transposes a matrix stored in a, ja, ia format.
c ---------------
c on entry:
c----------
c n	= dimension of A.
c job	= integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc (n,1,n+2,a,ja,ia,a,ja,ia(n+2)) 
c	  for any other normal usage, enter ipos=1.
c a	= real array of length nnz (nnz=number of nonzero elements in input 
c         matrix) containing the nonzero elements.
c ja	= integer array of length nnz containing the column positions
c 	  of the corresponding elements in a.
c ia	= integer of size n+1. ia(k) contains the position in a, ja of
c	  the beginning of the k-th row.
c
c on return:
c ---------- 
c output arguments:
c ao	= real array of size nzz containing the "a" part of the transpose
c jao	= integer array of size nnz containing the column indices.
c iao	= integer array of size n+1 containing the "ia" index array of
c	  the transpose. 
c
c----------------------------------------------------------------------- 
      call csrcsc2 (n,n,job,ipos,a,ja,ia,ao,jao,iao)
      end
c-----------------------------------------------------------------------
      subroutine csrcsc2 (n,n2,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n2+1),ja(*),jao(*)
      real*8  a(*),ao(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to      Compressed Sparse Column
c
c (transposition operation)   Not in place. 
c----------------------------------------------------------------------- 
c Rectangular version.  n is number of rows of CSR matrix,
c                       n2 (input) is number of columns of CSC matrix.
c----------------------------------------------------------------------- 
c -- not in place --
c this subroutine transposes a matrix stored in a, ja, ia format.
c ---------------
c on entry:
c----------
c n	= number of rows of CSR matrix.
c n2    = number of columns of CSC matrix.
c job	= integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc2 (n,n,1,n+2,a,ja,ia,a,ja,ia(n+2)) 
c	  for any other normal usage, enter ipos=1.
c a	= real array of length nnz (nnz=number of nonzero elements in input 
c         matrix) containing the nonzero elements.
c ja	= integer array of length nnz containing the column positions
c 	  of the corresponding elements in a.
c ia	= integer of size n+1. ia(k) contains the position in a, ja of
c	  the beginning of the k-th row.
c
c on return:
c ---------- 
c output arguments:
c ao	= real array of size nzz containing the "a" part of the transpose
c jao	= integer array of size nnz containing the column indices.
c iao	= integer array of size n+1 containing the "ia" index array of
c	  the transpose. 
c
c----------------------------------------------------------------------- 
c----------------- compute lengths of rows of transp(A) ----------------
      do 1 i=1,n2+1
         iao(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1 
            j = ja(k)+1
            iao(j) = iao(j)+1
 2       continue 
 3    continue
c---------- compute pointers from lengths ------------------------------
      iao(1) = ipos 
      do 4 i=1,n2
         iao(i+1) = iao(i) + iao(i+1)
 4    continue
c--------------- now do the actual copying ----------------------------- 
      do 6 i=1,n
         do 62 k=ia(i),ia(i+1)-1 
            j = ja(k) 
            next = iao(j)
            if (job .eq. 1)  ao(next) = a(k)
            jao(next) = i
            iao(j) = next+1
 62      continue
 6    continue
c-------------------------- reshift iao and leave ---------------------- 
      do 7 i=n2,1,-1
         iao(i+1) = iao(i)
 7    continue
      iao(1) = ipos
c--------------- end of csrcsc2 ---------------------------------------- 
c-----------------------------------------------------------------------
      end

c-----------------------------------------------------------------------
      subroutine atmux (n, x, y, a, ja, ia)
      real*8 x(*), y(*), a(*) 
      integer n, ia(*), ja(*)
c-----------------------------------------------------------------------
c         transp( A ) times a vector
c----------------------------------------------------------------------- 
c multiplies the transpose of a matrix by a vector when the original
c matrix is stored in compressed sparse row storage. Can also be
c viewed as the product of a matrix by a vector when the original
c matrix is stored in the compressed sparse column format.
c-----------------------------------------------------------------------
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=transp(A)*x
c
c-----------------------------------------------------------------------
c     local variables 
c
      integer i, k 
c-----------------------------------------------------------------------
c
c     zero out output vector
c 
      do 1 i=1,n
         y(i) = 0.0
 1    continue
c
c loop over the rows
c
      do 100 i = 1,n
         do 99 k=ia(i), ia(i+1)-1 
            y(ja(k)) = y(ja(k)) + x(i)*a(k)
 99      continue
 100  continue
c
      return
c-------------end-of-atmux---------------------------------------------- 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine atmuxr (m, n, x, y, a, ja, ia)
      real*8 x(*), y(*), a(*) 
      integer m, n, ia(*), ja(*)
c-----------------------------------------------------------------------
c         transp( A ) times a vector, A can be rectangular
c----------------------------------------------------------------------- 
c See also atmux.  The essential difference is how the solution vector
c is initially zeroed.  If using this to multiply rectangular CSC 
c matrices by a vector, m number of rows, n is number of columns.
c-----------------------------------------------------------------------
c
c on entry:
c----------
c m     = column dimension of A
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=transp(A)*x
c
c-----------------------------------------------------------------------
c     local variables 
c
      integer i, k 
c-----------------------------------------------------------------------
c
c     zero out output vector
c 
      do 1 i=1,m
         y(i) = 0.0
 1    continue
c
c loop over the rows
c
      do 100 i = 1,n
         do 99 k=ia(i), ia(i+1)-1 
            y(ja(k)) = y(ja(k)) + x(i)*a(k)
 99      continue
 100  continue
c
      return
c-------------end-of-atmuxr--------------------------------------------- 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine readmt_c (nmax,nzmax,job,fname,a,ja,ia,rhs,nrhs,
     *     guesol,nrow,ncol,nnz,title,key,type,ierr)
c-----------------------------------------------------------------------
c     this  subroutine reads  a boeing/harwell matrix,  given the
c     corresponding file. handles right hand sides in full format
c     only (no sparse right hand sides). Also the matrix must  be
c     in assembled forms.
c     It differs from readmt, in that  the name of the file needs
c     to be passed, and then the file is opened and closed within
c     this routine.
c     Author: Youcef Saad - Date: Oct 31, 1989
c     updated Jul 20, 1998 by Irene Moulitsas
c-----------------------------------------------------------------------
c on entry:
c---------
c nmax 	 =  max column dimension  allowed for matrix. The array ia should 
c	    be of length at least ncol+1 (see below) if job.gt.0
c nzmax	 = max number of nonzeros elements allowed. the arrays a, 
c          and ja should be of length equal to nnz (see below) if these
c          arrays are to be read (see job).
c          
c job	 = integer to indicate what is to be read. (note: job is an
c          input and output parameter, it can be modified on return)
c          job = 0    read the values of ncol, nrow, nnz, title, key,
c                     type and return. matrix is not read and arrays
c                     a, ja, ia, rhs are not touched.
c          job = 1    read srtucture only, i.e., the arrays ja and ia.
c          job = 2    read matrix including values, i.e., a, ja, ia
c          job = 3    read matrix and right hand sides: a,ja,ia,rhs.
c		      rhs may contain initial guesses and exact 
c                     solutions appended to the actual right hand sides.
c		      this will be indicated by the output parameter
c                     guesol [see below]. 
c                     
c fname = name of the file where to read the matrix from.
c
c nrhs   = integer. nrhs is an input as well as ouput parameter.
c          at input nrhs contains the total length of the array rhs.
c          See also ierr and nrhs in output parameters.
c
c on return:
c---------- 
c job    = on return job may be modified to the highest job it could
c          do: if job=2 on entry but no matrix values are available it
c          is reset to job=1 on return. Similarly of job=3 but no rhs 
c          is provided then it is rest to job=2 or job=1 depending on 
c          whether or not matrix values are provided.
c          Note that no error message is triggered (i.e. ierr = 0 
c          on return in these cases. It is therefore important to
c          compare the values of job on entry and return ).
c
c a	 = the a matrix in the a, ia, ja (column) storage format
c ja 	 = column number of element a(i,j) in array a.
c ia     = pointer  array. ia(i) points to the beginning of column i.
c
c rhs    = real array of size nrow + 1 if available (see job)
c
c nrhs   = integer containing the number of right-hand sides found
c          each right hand side may be accompanied with an intial guess
c          and also the exact solution.
c
c guesol = a 2-character string indicating whether an initial guess 
c          (1-st character) and / or the exact solution (2-nd
c          character) is provided with the right hand side.
c	   if the first character of guesol is 'G' it means that an
c          an intial guess is provided for each right-hand side.
c          These are appended to the right hand-sides in the array rhs.
c	   if the second character of guesol is 'X' it means that an
c          exact solution is provided for each right-hand side.
c          These are  appended to the right hand-sides 
c          and the initial guesses (if any) in the array rhs.
c
c nrow   = number of rows in matrix
c ncol	 = number of columns in matrix 
c nnz	 = number of nonzero elements in A. This info is returned
c          even if there is not enough space in a, ja, ia, in order
c          to determine the minimum storage needed. 
c
c title  = character*72 = title of matrix test ( character a*72). 
c key    = character*8  = key of matrix 
c type   = charatcer*3  = type of matrix.
c          for meaning of title, key and type refer to documentation 
c          Harwell/Boeing matrices.
c
c ierr   = integer used for error messages 
c         * ierr  =  0 means that  the matrix has been read normally. 
c         * ierr  =  1 means that  the array matrix could not be read 
c         because ncol+1 .gt. nmax
c         * ierr  =  2 means that  the array matrix could not be read 
c         because nnz .gt. nzmax 
c         * ierr  =  3 means that  the array matrix could not be read 
c         because both (ncol+1 .gt. nmax) and  (nnz .gt. nzmax )
c         * ierr  =  4 means that  the right hand side (s) initial 
c         guesse (s) and exact solution (s)   could  not be
c         read because they are stored in sparse format (not handled
c         by this routine ...) 
c         * ierr  =  5 means that the right-hand-sides, initial guesses
c         and exact solutions could not be read because the length of 
c         rhs as specified by the input value of nrhs is not 
c         sufficient to store them. The rest of the matrix may have
c         been read normally.
c 
c Notes:
c-------
c 1) This routine can be interfaced with the C language, since only
c    the name of the  file needs to be passed and no iounti number.
c
c 2) Refer to the  documentation on  the Harwell-Boeing formats for
c    details on the format assumed by readmt.
c    We summarize the format here for convenience.
c  
c    a) all  lines in  inout are  assumed to be 80  character long.
c    b) the file  consists of a header followed by the block of the
c       column start  pointers followed  by  the  block  of the row
c       indices,  followed  by  the  block  of the  real values and
c       finally the  numerical  values of  the right-hand-side if a
c       right hand side is supplied. 
c    c) the file starts by a header which contains four lines if no
c       right hand side is supplied and five lines otherwise.
c       * first  line  contains  the  title  (72  characters  long)
c         followed  by  the  8-character  identifier (name  of  the
c         matrix, called key) [ A72,A8 ]
c       * second line  contains the number of lines for each of the
c         following data blocks (4 of them) and the total number of
c         lines excluding the header.  [5i4]
c       * the   third  line  contains  a  three   character  string
c         identifying the type of  matrices as they  are referenced 
c         in  the Harwell  Boeing documentation [e.g., rua, rsa,..]
c         and the number of rows, columns, nonzero entries.
c         [A3,11X,4I14]
c       * The fourth  line contains the variable fortran format for
c         the following data blocks. [2A16,2A20] 
c       * The fifth  line is  only present if  right-hand-sides are
c         supplied. It  consists  of  three  one  character-strings
c         containing the  storage  format for the  right-hand-sides
c         ('F'= full,'M'=sparse=same as matrix), an  initial  guess
c         indicator  ('G' for yes),  an  exact  solution  indicator
c         ('X' for yes), followed by the number of right-hand-sides
c         and then the number of row indices.  [A3,11X,2I14] 
c     d) The three  following blocks follow the header as described
c        above.
c     e) In case the right hand-side are in sparse formats then the
c        fourth  block  uses the  same  storage  format as  for the 
c        matrix to  describe  the NRHS right  hand  sides provided,
c        with a column being replaced by a right hand side.
c-----------------------------------------------------------------------
      character title*72, key*8, type*3, ptrfmt*16, indfmt*16,
     &     valfmt*20, rhsfmt*20, rhstyp*3, guesol*2
      integer totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol,
     &     nnz, neltvl, nrhs, nmax, nzmax, nrwindx
      integer ia (nmax+1), ja (nzmax) 
      real*8 a(nzmax), rhs(*) 
      character fname*100
c-----------------------------------------------------------------------
      ierr = 0
      lenrhs = nrhs
c     
      iounit=15
      open(iounit,file = fname)
      read (iounit,10) title, key, totcrd, ptrcrd, indcrd, valcrd, 
     &     rhscrd, type, nrow, ncol, nnz, neltvl, ptrfmt, indfmt, 
     &     valfmt, rhsfmt
 10   format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
c     
      if (rhscrd .gt. 0) read (iounit,11) rhstyp, nrhs, nrwindx
 11   format (a3,11x,i14,i14)
c     
c     anything else to read ?
c     
      if (job .le. 0) goto 12

c     ---- check whether matrix is readable ------ 
      n = ncol
      if (ncol .gt. nmax) ierr = 1
      if (nnz .gt. nzmax) ierr = ierr + 2
      if (ierr .ne. 0) goto 12

c     ---- read pointer and row numbers ---------- 
      read (iounit,ptrfmt) (ia (i), i = 1, n+1)
      read (iounit,indfmt) (ja (i), i = 1, nnz)

c     --- reading values of matrix if required....
      if (job .le. 1)  goto 12
c     --- and if available ----------------------- 
      if (valcrd .le. 0) then
	 job = 1
	 goto 12
      endif
      read (iounit,valfmt) (a(i), i = 1, nnz)

c     --- reading rhs if required ---------------- 
      if (job .le. 2)  goto 12
c     --- and if available ----------------------- 
      if ( rhscrd .le. 0) then
	 job = 2
	 goto 12
      endif
c     
c     --- read right-hand-side.-------------------- 
c     
      if (rhstyp(1:1) .eq. 'M') then 
         ierr = 4
         goto 12
      endif
c     
      guesol = rhstyp(2:3) 
c     
      nvec = 1 
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') nvec=nvec+1
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') nvec=nvec+1
c     
      len = nrhs*nrow 
c     
      if (len*nvec .gt. lenrhs) then
        ierr = 5
        goto 12
      endif
c     
c     read right-hand-sides
c     
      next = 1
      iend = len
      read(iounit,rhsfmt) (rhs(i), i = next, iend)
c     
c     read initial guesses if available
c     
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') then
        next = next+len
        iend = iend+ len
        read(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c     
c     read exact solutions if available
c     
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') then
        next = next+len
        iend = iend+ len
        read(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c     
 12   close(iounit)
      return
c---------end-of-readmt_c-------------------------------------------------
c-------------------------------------------------------------------------
      end

c-------------------------------------------------------------------------
      subroutine csort (n,a,ja,ia,iwork,values) 
      logical values
      integer n, ja(*), ia(n+1), iwork(*) 
      real*8 a(*) 
c-----------------------------------------------------------------------
c This routine sorts the elements of  a matrix (stored in Compressed
c Sparse Row Format) in increasing order of their column indices within 
c each row. It uses a form of bucket sort with a cost of O(nnz) where
c nnz = number of nonzero elements. 
c requires an integer work array of length 2*nnz.  
c-----------------------------------------------------------------------
c on entry:
c--------- 
c n     = the row dimension of the matrix
c a     = the matrix A in compressed sparse row format.
c ja    = the array of column indices of the elements in array a.
c ia    = the array of pointers to the rows. 
c iwork = integer work array of length max ( n+1, 2*nnz ) 
c         where nnz = (ia(n+1)-ia(1))  ) .
c values= logical indicating whether or not the real values a(*) must 
c         also be permuted. if (.not. values) then the array a is not
c         touched by csort and can be a dummy array. 
c 
c on return:
c----------
c the matrix stored in the structure a, ja, ia is permuted in such a
c way that the column indices are in increasing order within each row.
c iwork(1:nnz) contains the permutation used  to rearrange the elements.
c----------------------------------------------------------------------- 
c Y. Saad - Feb. 1, 1991.
c revised by Zhongze Li - June 13th, 2001
c-----------------------------------------------------------------------
c local variables
      integer i, k, j, ifirst, offset, nnz, next  
c
c count the number of elements in each column
c
      offset = 1 - ia(1)
      do 1 i=1,n+1
         iwork(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1 
            j = ja(k+offset)+1+offset
            iwork(j) = iwork(j)+1
 2       continue 
 3    continue
c
c compute pointers from lengths. 
c
      iwork(1) = 1
      do 4 i=1,n
         iwork(i+1) = iwork(i) + iwork(i+1)
 4    continue
c 
c get the positions of the nonzero elements in order of columns.
c
      ifirst = ia(1); 
      nnz = ia(n+1)-ifirst
      do 5 i=1,n
         do 51 k=ia(i),ia(i+1)-1 
            j = ja(k+offset)+offset 
            next = iwork(j) 
            iwork(nnz+next) = k+offset
            iwork(j) = next+1
 51      continue
 5    continue
c
c convert to coordinate format
c 
      do 6 i=1, n
         do 61 k=ia(i), ia(i+1)-1 
            iwork(k+offset) = i
 61      continue
 6    continue
c
c loop to find permutation: for each element find the correct 
c position in (sorted) arrays a, ja. Record this in iwork. 
c 
      do 7 k=1, nnz
         ko = iwork(nnz+k) 
         irow = iwork(ko)
         next = ia(irow)
c
c the current element should go in next position in row. iwork
c records this position. 
c 
         iwork(ko) = next+offset
         ia(irow)  = next+1
 7       continue
c
c perform an in-place permutation of the  arrays.
c 
         call ivperm (nnz, ja, iwork) 
         if (values) call dvperm (nnz, a, iwork) 
c
c reshift the pointers of the original matrix back.
c 
      do 8 i=n,1,-1
         ia(i+1) = ia(i)
 8    continue
      ia(1) = ifirst 
c
      return 
c---------------end-of-csort-------------------------------------------- 
c-----------------------------------------------------------------------
      end

c----------------------------------------------------------------------- 
      subroutine dvperm (n, x, perm) 
      integer n, perm(n) 
      real*8 x(n)
c-----------------------------------------------------------------------
c this subroutine performs an in-place permutation of a real vector x 
c according to the permutation array perm(*), i.e., on return, 
c the vector x satisfies,
c
c	x(perm(j)) :== x(j), j=1,2,.., n
c
c-----------------------------------------------------------------------
c on entry:
c---------
c n 	= length of vector x.
c perm 	= integer array of length n containing the permutation  array.
c x	= input vector
c
c on return:
c---------- 
c x	= vector x permuted according to x(perm(*)) :=  x(*)
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
c local variables 
      real*8 tmp, tmp1
c
      init      = 1
      tmp	= x(init)	
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
c     
c loop
c 
 6    k = k+1
c
c save the chased element --
c 
      tmp1	  = x(ii) 
      x(ii)     = tmp
      next	  = perm(ii) 
      if (next .lt. 0 ) goto 65
c     
c test for end 
c
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next 
c
c end loop 
c
      goto 6
c
c reinitilaize cycle --
c
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp	= x(init)
      ii	= perm(init)
      perm(init)=-perm(init)
      goto 6
c     
 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue 
c     
      return
c-------------------end-of-dvperm--------------------------------------- 
c-----------------------------------------------------------------------
      end

c-----------------------------------------------------------------------
      subroutine ivperm (n, ix, perm) 
      integer n, perm(n), ix(n)
c-----------------------------------------------------------------------
c this subroutine performs an in-place permutation of an integer vector 
c ix according to the permutation array perm(*), i.e., on return, 
c the vector x satisfies,
c
c	ix(perm(j)) :== ix(j), j=1,2,.., n
c
c-----------------------------------------------------------------------
c on entry:
c---------
c n 	= length of vector x.
c perm 	= integer array of length n containing the permutation  array.
c ix	= input vector
c
c on return:
c---------- 
c ix	= vector x permuted according to ix(perm(*)) :=  ix(*)
c
c----------------------------------------------------------------------c
c           Y. Saad, Sep. 21 1989                                      c
c----------------------------------------------------------------------c
c local variables
      integer tmp, tmp1
c
      init      = 1
      tmp	= ix(init)	
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
c     
c loop
c 
 6    k = k+1
c
c save the chased element --
c 
      tmp1	  = ix(ii) 
      ix(ii)     = tmp
      next	  = perm(ii) 
      if (next .lt. 0 ) goto 65
c     
c test for end 
c
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next 
c
c end loop 
c
      goto 6
c
c reinitilaize cycle --
c
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp	= ix(init)
      ii	= perm(init)
      perm(init)=-perm(init)
      goto 6
c     
 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue 
c     
      return
c-------------------end-of-ivperm--------------------------------------- 
c-----------------------------------------------------------------------
      end


c----------------------------------------------------------------------- 
      subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),job
      real*8 a(*),ao(*) 
c-----------------------------------------------------------------------
c this subroutine permutes the rows of a matrix in CSR format. 
c rperm  computes B = P A  where P is a permutation matrix.  
c the permutation P is defined through the array perm: for each j, 
c perm(j) represents the destination row number of row number j. 
c Youcef Saad -- recoded Jan 28, 1991.
c-----------------------------------------------------------------------
c on entry:
c----------
c n 	= dimension of the matrix
c a, ja, ia = input matrix in csr format
c perm 	= integer array of length nrow containing the permutation arrays
c	  for the rows: perm(i) is the destination of row i in the
c         permuted matrix. 
c         ---> a(i,j) in the original matrix becomes a(perm(i),j) 
c         in the output  matrix.
c
c job	= integer indicating the work to be done:
c 		job = 1	permute a, ja, ia into ao, jao, iao 
c                       (including the copying of real values ao and
c                       the array iao).
c 		job .ne. 1 :  ignore real values.
c                     (in which case arrays a and ao are not needed nor
c                      used).
c
c------------
c on return: 
c------------ 
c ao, jao, iao = input matrix in a, ja, ia format
c note : 
c        if (job.ne.1)  then the arrays a and ao are not used.
c----------------------------------------------------------------------c
c           Y. Saad, May  2, 1990                                      c
c----------------------------------------------------------------------c
      logical values
      values = (job .eq. 1) 
c     
c     determine pointers for output matix. 
c     
      do 50 j=1,nrow
         i = perm(j)
         iao(i+1) = ia(j+1) - ia(j)
 50   continue
c
c get pointers from lengths
c
      iao(1) = 1
      do 51 j=1,nrow
         iao(j+1)=iao(j+1)+iao(j)
 51   continue
c
c copying 
c
      do 100 ii=1,nrow
c
c old row = ii  -- new row = iperm(ii) -- ko = new pointer
c        
         ko = iao(perm(ii)) 
         do 60 k=ia(ii), ia(ii+1)-1 
            jao(ko) = ja(k) 
            if (values) ao(ko) = a(k)
            ko = ko+1
 60      continue
 100  continue
c
      return
c---------end-of-rperm ------------------------------------------------- 
c-----------------------------------------------------------------------
      end

      subroutine cperm (nrow,a,ja,ia,ao,jao,iao,perm,job) 
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*), job
      real*8 a(*), ao(*) 
c-----------------------------------------------------------------------
c this subroutine permutes the columns of a matrix a, ja, ia.
c the result is written in the output matrix  ao, jao, iao.
c cperm computes B = A P, where  P is a permutation matrix
c that maps column j into column perm(j), i.e., on return 
c      a(i,j) becomes a(i,perm(j)) in new matrix 
c Y. Saad, May 2, 1990 / modified Jan. 28, 1991. 
c-----------------------------------------------------------------------
c on entry:
c----------
c nrow 	= row dimension of the matrix
c
c a, ja, ia = input matrix in csr format. 
c
c perm	= integer array of length ncol (number of columns of A
c         containing the permutation array  the columns: 
c         a(i,j) in the original matrix becomes a(i,perm(j))
c         in the output matrix.
c
c job	= integer indicating the work to be done:
c 		job = 1	permute a, ja, ia into ao, jao, iao 
c                       (including the copying of real values ao and
c                       the array iao).
c 		job .ne. 1 :  ignore real values ao and ignore iao.
c
c------------
c on return: 
c------------ 
c ao, jao, iao = input matrix in a, ja, ia format (array ao not needed)
c
c Notes:
c------- 
c 1. if job=1 then ao, iao are not used.
c 2. This routine is in place: ja, jao can be the same. 
c 3. If the matrix is initially sorted (by increasing column number) 
c    then ao,jao,iao  may not be on return. 
c 
c----------------------------------------------------------------------c
c local parameters:
      integer k, i, nnz
c
      nnz = ia(nrow+1)-1
      do 100 k=1,nnz
         jao(k) = perm(ja(k)) 
 100  continue
c
c     done with ja array. return if no need to touch values.
c
      if (job .ne. 1) return
c
c else get new pointers -- and copy values too.
c 
      do 1 i=1, nrow+1
         iao(i) = ia(i)
 1    continue
c
      do 2 k=1, nnz
         ao(k) = a(k)
 2    continue
c
      return
c---------end-of-cperm-------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c
      subroutine wreadmtc (nmax,nzmax,job,fname,length,a,ja,ia,rhs,nrhs,
     *     guesol,nrow,ncol,nnz,title,key,type,ierr)
c-----------------------------------------------------------------------
c     this  subroutine reads  a boeing/harwell matrix,  given the
c     corresponding file. handles right hand sides in full format
c     only (no sparse right hand sides). Also the matrix must  be
c     in assembled forms.
c     It differs from readmt, in that  the name of the file needs
c     to be passed, and then the file is opened and closed within
c     this routine.
c     Author: Youcef Saad - Date: Oct 31, 1989
c     updated Feb 6th, 2001 by Zhongze Li
c-----------------------------------------------------------------------
      character fname*100, title*72, key*8, type*3, guesol*2, fname1*100
      integer nrow, ncol, nnz, nrhs, nmax, nzmax, job, length, ierr
      integer ia (nmax+1), ja (nzmax)
      real*8 a(nzmax), rhs(*)

      fname1 = " "
      fname1(1:length) = fname(1:length)
c      print *,nmax,nzmax,job,fname,length,nrhs,
c     *     nrow,ncol,nnz,title,key,type,ierr
      call readmtc(nmax,nzmax,job,fname1,a,ja,ia,rhs,nrhs,
     *     guesol,nrow,ncol,nnz,title,key,type,ierr)
c      print *,nmax,nzmax,job,fname,length,nrhs,
c     *     nrow,ncol,nnz,title,key,type,ierr
      return
      end

      subroutine readmtc (nmax,nzmax,job,fname,a,ja,ia,rhs,nrhs,
     *     guesol,nrow,ncol,nnz,title,key,type,ierr)
c-----------------------------------------------------------------------
c     this  subroutine reads  a boeing/harwell matrix,  given the
c     corresponding file. handles right hand sides in full format
c     only (no sparse right hand sides). Also the matrix must  be
c     in assembled forms.
c     It differs from readmt, in that  the name of the file needs
c     to be passed, and then the file is opened and closed within
c     this routine.
c     Author: Youcef Saad - Date: Oct 31, 1989
c     updated Jul 20, 1998 by Irene Moulitsas
c-----------------------------------------------------------------------
c on entry:
c---------
c nmax   =  max column dimension  allowed for matrix. The array ia should
c           be of length at least ncol+1 (see below) if job.gt.0
c nzmax  = max number of nonzeros elements allowed. the arrays a,
c          and ja should be of length equal to nnz (see below) if these
c          arrays are to be read (see job).
c
c job    = integer to indicate what is to be read. (note: job is an
c          input and output parameter, it can be modified on return)
c          job = 0    read the values of ncol, nrow, nnz, title, key,
c                     type and return. matrix is not read and arrays
c                     a, ja, ia, rhs are not touched.
c          job = 1    read srtucture only, i.e., the arrays ja and ia.
c          job = 2    read matrix including values, i.e., a, ja, ia
c          job = 3    read matrix and right hand sides: a,ja,ia,rhs.
c                     rhs may contain initial guesses and exact
c                     solutions appended to the actual right hand sides.
c                     this will be indicated by the output parameter
c                     guesol [see below].
c
c fname = name of the file where to read the matrix from.
c
c nrhs   = integer. nrhs is an input as well as ouput parameter.
c          at input nrhs contains the total length of the array rhs.
c          See also ierr and nrhs in output parameters.
c
c on return:
c----------
c job    = on return job may be modified to the highest job it could
c          do: if job=2 on entry but no matrix values are available it
c          is reset to job=1 on return. Similarly of job=3 but no rhs
c          is provided then it is rest to job=2 or job=1 depending on
c          whether or not matrix values are provided.
c          Note that no error message is triggered (i.e. ierr = 0
c          on return in these cases. It is therefore important to
c          compare the values of job on entry and return ).
c
c a      = the a matrix in the a, ia, ja (column) storage format
c ja     = column number of element a(i,j) in array a.
c ia     = pointer  array. ia(i) points to the beginning of column i.
c
c rhs    = real array of size nrow + 1 if available (see job)
c
c nrhs   = integer containing the number of right-hand sides found
c          each right hand side may be accompanied with an intial guess
c          and also the exact solution.
c
c guesol = a 2-character string indicating whether an initial guess
c          (1-st character) and / or the exact solution (2-nd
c          character) is provided with the right hand side.
c          if the first character of guesol is 'G' it means that an
c          an intial guess is provided for each right-hand side.
c          These are appended to the right hand-sides in the array rhs.
c          if the second character of guesol is 'X' it means that an
c          exact solution is provided for each right-hand side.
c          These are  appended to the right hand-sides
c          and the initial guesses (if any) in the array rhs.
c
c nrow   = number of rows in matrix
c ncol   = number of columns in matrix
c nnz    = number of nonzero elements in A. This info is returned
c          even if there is not enough space in a, ja, ia, in order
c          to determine the minimum storage needed.
c
c title  = character*72 = title of matrix test ( character a*72).
c key    = character*8  = key of matrix
c type   = charatcer*3  = type of matrix.
c          for meaning of title, key and type refer to documentation
c          Harwell/Boeing matrices.
c
c ierr   = integer used for error messages
c         * ierr  =  0 means that  the matrix has been read normally.
c         * ierr  =  1 means that  the array matrix could not be read
c         because ncol+1 .gt. nmax
c         * ierr  =  2 means that  the array matrix could not be read
c         because nnz .gt. nzmax
c         * ierr  =  3 means that  the array matrix could not be read
c         because both (ncol+1 .gt. nmax) and  (nnz .gt. nzmax )
c         * ierr  =  4 means that  the right hand side (s) initial
c         guesse (s) and exact solution (s)   could  not be
c         read because they are stored in sparse format (not handled
c         by this routine ...)
c         * ierr  =  5 means that the right-hand-sides, initial guesses
c         and exact solutions could not be read because the length of
c         rhs as specified by the input value of nrhs is not
c         sufficient to store them. The rest of the matrix may have
c         been read normally.
c
c Notes:
c-------
c 1) This routine can be interfaced with the C language, since only
c    the name of the  file needs to be passed and no iounti number.
c
c 2) Refer to the  documentation on  the Harwell-Boeing formats for
c    details on the format assumed by readmt.
c    We summarize the format here for convenience.
c
c    a) all  lines in  inout are  assumed to be 80  character long.
c    b) the file  consists of a header followed by the block of the
c       column start  pointers followed  by  the  block  of the row
c       indices,  followed  by  the  block  of the  real values and
c       finally the  numerical  values of  the right-hand-side if a
c       right hand side is supplied.
c    c) the file starts by a header which contains four lines if no
c       right hand side is supplied and five lines otherwise.
c       * first  line  contains  the  title  (72  characters  long)
c         followed  by  the  8-character  identifier (name  of  the
c         matrix, called key) [ A72,A8 ]
c       * second line  contains the number of lines for each of the
c         following data blocks (4 of them) and the total number of
c         lines excluding the header.  [5i4]
c       * the   third  line  contains  a  three   character  string
c         identifying the type of  matrices as they  are referenced
c         in  the Harwell  Boeing documentation [e.g., rua, rsa,..]
c         and the number of rows, columns, nonzero entries.
c         [A3,11X,4I14]
c       * The fourth  line contains the variable fortran format for
c         the following data blocks. [2A16,2A20]
c       * The fifth  line is  only present if  right-hand-sides are
c         supplied. It  consists  of  three  one  character-strings
c         containing the  storage  format for the  right-hand-sides
c         ('F'= full,'M'=sparse=same as matrix), an  initial  guess
c         indicator  ('G' for yes),  an  exact  solution  indicator
c         ('X' for yes), followed by the number of right-hand-sides
c         and then the number of row indices.  [A3,11X,2I14]
c     d) The three  following blocks follow the header as described
c        above.
c     e) In case the right hand-side are in sparse formats then the
c        fourth  block  uses the  same  storage  format as  for the
c        matrix to  describe  the NRHS right  hand  sides provided,
c        with a column being replaced by a right hand side.
c-----------------------------------------------------------------------
      character title*72, key*8, type*3, ptrfmt*16, indfmt*16,
     &     valfmt*20, rhsfmt*20, rhstyp*3, guesol*2
      integer totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol,
     &     nnz, neltvl, nrhs, nmax, nzmax, nrwindx
      integer ia (nmax+1), ja (nzmax)
      real*8 a(nzmax), rhs(*)
      character fname*100
c-----------------------------------------------------------------------
      ierr = 0
      lenrhs = nrhs
c
      iounit=15
      open(iounit,file = fname)
      read (iounit,10) title, key, totcrd, ptrcrd, indcrd, valcrd,
     &     rhscrd, type, nrow, ncol, nnz, neltvl, ptrfmt, indfmt,
     &     valfmt, rhsfmt
 10   format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
c
      if (rhscrd .gt. 0) read (iounit,11) rhstyp, nrhs, nrwindx
 11   format (a3,11x,i14,i14)
c
c     anything else to read ?
c
      if (job .le. 0) goto 12
c     ---- check whether matrix is readable ------
      n = ncol
      if (ncol .gt. nmax) ierr = 1
      if (nnz .gt. nzmax) ierr = ierr + 2
      if (ierr .ne. 0) goto 12
c     ---- read pointer and row numbers ----------
      read (iounit,ptrfmt) (ia (i), i = 1, n+1)
      read (iounit,indfmt) (ja (i), i = 1, nnz)
c     --- reading values of matrix if required....
      if (job .le. 1)  goto 12
c     --- and if available -----------------------
      if (valcrd .le. 0) then
         job = 1
         goto 12
      endif
      read (iounit,valfmt) (a(i), i = 1, nnz)
c     --- reading rhs if required ----------------
      if (job .le. 2)  goto 12
c     --- and if available -----------------------
      if ( rhscrd .le. 0) then
         job = 2
         goto 12
      endif
c
c     --- read right-hand-side.--------------------
c
      if (rhstyp(1:1) .eq. 'M') then
         ierr = 4
         goto 12
      endif
c
      guesol = rhstyp(2:3)
c
      nvec = 1
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') nvec=nvec+1
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') nvec=nvec+1
c
      len = nrhs*nrow
c
      if (len*nvec .gt. lenrhs) then
        ierr = 5
        goto 12
      endif
c
c     read right-hand-sides
c
      next = 1
      iend = len
      read(iounit,rhsfmt) (rhs(i), i = next, iend)
c
c     read initial guesses if available
c
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') then
        next = next+len
        iend = iend+ len
        read(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c
c     read exact solutions if available
c
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') then
        next = next+len
        iend = iend+ len
        read(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c
 12   close(iounit)
      return
c---------end-of-readmt_c-------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
c     some routines extracted from SPARSKIT2 and BLAS.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      function distdot(n,x,ix,y,iy)
      integer n, ix, iy
      real*8 distdot, x(*), y(*), ddot
      external ddot
      distdot = ddot(n,x,ix,y,iy)
      return
      end
c-----end-of-distdot
c-----------------------------------------------------------------------

c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c               REORDERING ROUTINES -- LEVEL SET BASED ROUTINES        c
c----------------------------------------------------------------------c
c dblstr   : doubled stripe partitioner 
c rdis     : recursive dissection partitioner
c dse2way  : distributed site expansion usuing sites from dblstr 
c dse      : distributed site expansion usuing sites from rdis
c------------- utility routines ----------------------------------------- 
c BFS      : Breadth-First search traversal algorithm 
c add_lvst : routine to add a level -- used by BFS 
c stripes  : finds the level set structure
c stripes0 : finds a trivial one-way partitioning from level-sets 
c perphn   : finds a pseudo-peripheral node and performs a BFS from it.
c mapper4  : routine used by dse and dse2way to do center expansion
c get_domns: routine to find subdomaine from linked lists found by 
c            mapper4. 
c add_lk   : routine to add entry to linked list -- used by mapper4. 
c find_ctr : routine to locate an approximate center of a subgraph. 
c rversp   : routine to reverse a given permutation (e.g., for RCMK)
c maskdeg  : integer function to compute the `masked' of a node
c-----------------------------------------------------------------------
      subroutine dblstr(n,ja,ia,ip1,ip2,nfirst,riord,ndom,map,mapptr,
     *     mask,levels,iwk) 
      implicit none
      integer ndom,ja(*),ia(*),ip1,ip2,nfirst,riord(*),map(*),mapptr(*),
     *     mask(*),levels(*),iwk(*),nextdom
c-----------------------------------------------------------------------
c     this routine does a two-way partitioning of a graph using 
c     level sets recursively. First a coarse set is found by a
c     simple cuthill-mc Kee type algorithm. Them each of the large
c     domains is further partitioned into subsets using the same 
c     technique. The ip1 and ip2 parameters indicate the desired number 
c     number of partitions 'in each direction'. So the total number of
c     partitions on return ought to be equal (or close) to ip1*ip2 
c----------------------parameters----------------------------------------
c on entry: 
c---------
c n      = row dimension of matrix == number of vertices in graph
c ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data
c          structure)
c ip1    = integer indicating the number of large partitions ('number of
c          paritions in first direction') 
c ip2    = integer indicating the number of smaller partitions, per 
c          large partition, ('number of partitions in second direction') 
c nfirst = number of nodes in the first level that is input in riord 
c riord  = (also an ouput argument). on entry riord contains the labels  
c          of the nfirst nodes that constitute the first level.   
c on return:
c-----------
c ndom   = total number of partitions found 
c map    = list of nodes listed partition by partition from partition 1
c          to paritition ndom.
c mapptr = pointer array for map. All nodes from position 
c          k1=mapptr(idom),to position k2=mapptr(idom+1)-1 in map belong
c          to partition idom.
c work arrays:
c-------------
c mask   = array of length n, used to hold the partition number of each 
c          node for the first (large) partitioning. 
c          mask is also used as a marker of  visited nodes. 
c levels = integer array of length .le. n used to hold the pointer 
c          arrays for the various level structures obtained from BFS. 
c 
c-----------------------------------------------------------------------
      integer n, j,idom,kdom,jdom,maskval,k,nlev,init,ndp1,numnod
      maskval = 1 
      do j=1, n
         mask(j) = maskval 
      enddo
      iwk(1) = 0 
      call BFS(n,ja,ia,nfirst,iwk,mask,maskval,riord,levels,nlev)      
c
c     init = riord(1) 
c     call perphn (ja,ia,mask,maskval,init,nlev,riord,levels) 
      call stripes (nlev,riord,levels,ip1,map,mapptr,ndom)
c-----------------------------------------------------------------------
      if (ip2 .eq. 1) return      
      ndp1 = ndom+1
c     
c     pack info into array iwk 
c 
      do j = 1, ndom+1
         iwk(j) = ndp1+mapptr(j)  
      enddo
      do j=1, mapptr(ndom+1)-1
         iwk(ndp1+j) = map(j) 
      enddo
      do idom=1, ndom 
         j = iwk(idom) 
         numnod = iwk(idom+1) - iwk(idom) 
         init = iwk(j) 
         do k=j, iwk(idom+1)-1 
         enddo
      enddo

      do idom=1, ndom 
         do k=mapptr(idom),mapptr(idom+1)-1 
            mask(map(k)) = idom
         enddo
      enddo
      nextdom = 1 
c
c     jdom = counter for total number of (small) subdomains 
c 
      jdom = 1
      mapptr(jdom) = 1 
c----------------------------------------------------------------------- 
      do idom =1, ndom
         maskval = idom
         nfirst = 1
         numnod = iwk(idom+1) - iwk(idom) 
         j = iwk(idom) 
         init = iwk(j) 
         nextdom = mapptr(jdom) 
c  note:    old version uses iperm array 
         call perphn(numnod,ja,ia,init,mask,maskval,
     *        nlev,riord,levels)
c          
         call stripes (nlev,riord,levels,ip2,map(nextdom),
     *        mapptr(jdom),kdom)
c          
         mapptr(jdom) = nextdom
         do j = jdom,jdom+kdom-1
            mapptr(j+1) = nextdom + mapptr(j+1)-1
         enddo
         jdom = jdom + kdom
      enddo
c
      ndom = jdom - 1
      return
      end 
c-----------------------------------------------------------------------
      subroutine rdis(n,ja,ia,ndom,map,mapptr,mask,levels,size,iptr) 
      implicit none
c----- 
      integer n,ja(*),ia(*),ndom,map(*),mapptr(*),mask(*),levels(*),
     *     size(ndom+1),iptr(ndom+1)
c-----------------------------------------------------------------------
c     recursive dissection algorithm for partitioning.
c     initial graph is cut in two - then each time, the largest set
c     is cut in two until we reach desired number of domains.
c-----------------------------------------------------------------------
c     input
c     n, ja, ia = graph 
c     ndom      = desired number of subgraphs
c     output 
c     ------
c     map, mapptr  = pointer array data structure for domains. 
c     if k1 = mapptr(i), k2=mapptr(i+1)-1 then 
c     map(k1:k2) = points in domain number i
c     work arrays:
c    -------------
c    mask(1:n)    integer 
c    levels(1:n)  integer 
c    size(1:ndom) integer 
c    iptr(1:ndom) integer 
c----------------------------------------------------------------------- 
      integer idom,maskval,k,nlev,init,nextsiz,wantsiz,lev,ko,
     *     maxsiz,j,nextdom  
c-----------------------------------------------------------------------
      idom = 1
c-----------------------------------------------------------------------
c     size(i) = size of domnain  i
c     iptr(i)  = index of first element of domain i
c-----------------------------------------------------------------------
      size(idom) = n
      iptr(idom) = 1
      do j=1, n
         mask(j) = 1 
      enddo
c     
c     domain loop
c
 1    continue
c
c     select domain with largest size
c     
      maxsiz = 0 
      do j=1, idom
         if (size(j) .gt. maxsiz) then
            maxsiz = size(j)
            nextdom = j
         endif
      enddo
c
c     do a Prphn/ BFS on nextdom
c     
      maskval = nextdom
      init = iptr(nextdom) 
      call perphn(n,ja,ia,init,mask,maskval,nlev,map,levels) 
c
c     determine next subdomain
c
      nextsiz = 0
      wantsiz = maxsiz/2 
      idom = idom+1
      lev = nlev 
      do while (nextsiz .lt. wantsiz) 
         do k = levels(lev), levels(lev+1)-1
            mask(map(k)) = idom
        enddo
        nextsiz = nextsiz + levels(lev+1) - levels(lev) 
        lev = lev-1
      enddo
c
      size(nextdom) = size(nextdom) - nextsiz
      size(idom) = nextsiz
c
c     new initial point = last point of previous domain
c
       iptr(idom) = map(levels(nlev+1)-1) 
c       iptr(idom) = map(levels(lev)+1) 
c      iptr(idom) = 1 
c
c alternative 
c      lev = 1 
c      do while (nextsiz .lt. wantsiz) 
c         do k = levels(lev), levels(lev+1)-1
c            mask(map(k)) = idom
c         enddo
c         nextsiz = nextsiz + levels(lev+1) - levels(lev) 
c         lev = lev+1
c      enddo
c
c     set size of new domain and adjust previous one
c
c      size(idom) = nextsiz 
c      size(nextdom) = size(nextdom) - nextsiz 
c      iptr(idom) = iptr(nextdom) 
c      iptr(nextdom) = map(levels(lev))

      if (idom .lt. ndom) goto 1
c
c     domains found -- build data structure 
c     
      mapptr(1) = 1
      do idom=1, ndom 
         mapptr(idom+1) = mapptr(idom) + size(idom) 
      enddo
      do k=1, n
         idom = mask(k) 
         ko = mapptr(idom) 
         map(ko) = k
         mapptr(idom) = ko+1
      enddo
c
c     reset pointers
c     
      do j = ndom,1,-1
         mapptr(j+1) = mapptr(j) 
      enddo
      mapptr(1) = 1 
c
      return
      end 
c-----------------------------------------------------------------------
      subroutine dse2way(n,ja,ia,ip1,ip2,nfirst,riord,ndom,dom,idom,
     *     mask,jwk,link) 
c-----------------------------------------------------------------------
c     uses centers obtained from dblstr partition to get new partition
c----------------------------------------------------------------------- 
c     input: n, ja, ia   = matrix
c     nfirst = number of first points 
c     riord  = riord(1:nfirst) initial points 
c     output 
c     ndom   = number of domains
c     dom, idom = pointer array structure for domains. 
c     mask , jwk, link = work arrays,
c-----------------------------------------------------------------------
      implicit none 
      integer n, ja(*), ia(*), ip1, ip2, nfirst, riord(*), dom(*),
     *     idom(*), mask(*), jwk(n+1),ndom,link(*)  
c
c-----------------------------------------------------------------------
c     local variables
      integer i, mid,nsiz, maskval,init, outer, nouter, k
      call dblstr(n,ja,ia,ip1,ip2,nfirst,riord,ndom,dom,idom,mask,
     *     link,jwk)
c
      nouter = 3
c----------------------------------------------------------------------- 

      do outer =1, nouter 
c
c     set masks 
c
      do i=1, ndom
         do k=idom(i),idom(i+1)-1
            mask(dom(k)) = i
         enddo
      enddo
c
c     get centers 
c 
      do i =1, ndom
         nsiz = idom(i+1) - idom(i) 
         init = dom(idom(i))
         maskval = i 
c
c         use link for local riord -- jwk for other arrays -- 
c 
         call find_ctr(n,nsiz,ja,ia,init,mask,maskval,link, 
     *    jwk,mid,jwk(nsiz+1)) 
         riord(i) = mid 
      enddo
c
c     do level-set expansion from centers -- save previous diameter 
c 
      call mapper4(n,ja,ia,ndom,riord,jwk,mask,link) 
      call get_domns2(ndom,riord,link,jwk,dom,idom)
c----------------------------------------------------------------------- 
      enddo 
      return 
      end 
c-----------------------------------------------------------------------
      subroutine dse(n,ja,ia,ndom,riord,dom,idom,mask,jwk,link) 
      implicit none 
      integer n, ja(*), ia(n+1), ndom, riord(n), dom(*),
     *     idom(ndom+1),mask(n),jwk(2*ndom),link(n)  
c-----------------------------------------------------------------------
c     uses centers produced from rdis to get a new partitioning -- 
c     see calling sequence in rdis.. 
c     on entry:
c     n, ja, ia = graph 
c     ndom = number of desired subdomains
c     on return
c     dom, idom = 
c          pointer array structure for the ndom domains. 
c
c-----------------------------------------------------------------------
c     size of array jwk = 2*ndom
c     size array link
c     dom = array of size at least n
c     mask = array of size n.
c     link = same size as riord 
c     riord = n list of nodes listed in sequence for all subdomains 
c-----------------------------------------------------------------------
c     local variables
      integer i, mid, nsiz, maskval,init, outer, nouter, k 
c-----------------------------------------------------------------------
      nouter = 3
c 
      call rdis(n,ja,ia,ndom,dom,idom,mask,link,jwk,jwk(ndom+1)) 
c
c     initial points = 
c
      do outer =1, nouter 
c
c     set masks 
c
      do i=1, ndom
         do k=idom(i),idom(i+1)-1
            mask(dom(k)) = i
         enddo
      enddo
c
c     get centers 
c 
      do i =1, ndom
         nsiz = idom(i+1) - idom(i) 
         init = dom(idom(i))
         maskval = i 
c
c         use link for local riord -- jwk for other arrays -- 
c 
         call find_ctr(n,nsiz,ja,ia,init,mask,maskval,link, 
     *    jwk,mid,jwk(nsiz+1)) 
         riord(i) = mid 
      enddo
c
c     do level-set expansion from centers -- save previous diameter 
c 
      call mapper4(n,ja,ia,ndom,riord,jwk,mask,link) 
      call get_domns2(ndom,riord,link,jwk,dom,idom)
c----------------------------------------------------------------------- 
      enddo 
      return 
      end 
c----------------------------------------------------------------------- 
      subroutine BFS(n,ja,ia,nfirst,iperm,mask,maskval,riord,levels,
     *     nlev)
      implicit none 
      integer n,ja(*),ia(*),nfirst,iperm(n),mask(n+1),riord(*),
     *     levels(*),
     *     nlev,maskval 
c-----------------------------------------------------------------------
c finds the level-structure (breadth-first-search or CMK) ordering for a
c given sparse matrix. Uses add_lvst. Allows an set of nodes to be 
c the initial level (instead of just one node). 
c-------------------------parameters------------------------------------
c on entry:
c---------
c     n      = number of nodes in the graph 
c     ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data
c     structure)
c     nfirst = number of nodes in the first level that is input in riord
c     iperm  = integer array indicating in which order to  traverse the graph
c     in order to generate all connected components. 
c     if iperm(1) .eq. 0 on entry then BFS will traverse the nodes
c     in the  order 1,2,...,n.
c     
c     riord  = (also an ouput argument). On entry riord contains the labels  
c     of the nfirst nodes that constitute the first level.      
c     
c     mask   = array used to indicate whether or not a node should be 
c     condidered in the graph. see maskval.
c     mask is also used as a marker of  visited nodes. 
c     
c     maskval= consider node i only when:  mask(i) .eq. maskval 
c     maskval must be .gt. 0. 
c     thus, to consider all nodes, take mask(1:n) = 1. 
c     maskval=1 (for example) 
c     
c     on return
c     ---------
c     mask   = on return mask is restored to its initial state. 
c     riord  = `reverse permutation array'. Contains the labels of the nodes
c     constituting all the levels found, from the first level to
c     the last. 
c     levels = pointer array for the level structure. If lev is a level
c     number, and k1=levels(lev),k2=levels(lev+1)-1, then
c     all the nodes of level number lev are:
c     riord(k1),riord(k1+1),...,riord(k2) 
c     nlev   = number of levels found
c-----------------------------------------------------------------------
c     
      integer j, ii, nod, istart, iend 
      logical permut
      permut = (iperm(1) .ne. 0) 
c     
c     start pointer structure to levels 
c     
      nlev   = 0 
c     
c     previous end
c     
      istart = 0 
      ii = 0
c     
c     current end 
c     
      iend = nfirst
c     
c     intialize masks to zero -- except nodes of first level -- 
c     
      do 12 j=1, nfirst 
c         print *,'mask',riord(j),j,nfirst
         mask(riord(j)) = 0 
 12   continue
c-----------------------------------------------------------------------
 13   continue 
c     
 1    nlev = nlev+1
      levels(nlev) = istart + 1
      call add_lvst (istart,iend,nlev,riord,ja,ia,mask,maskval) 
      if (istart .lt. iend) goto 1
 2    ii = ii+1 
      if (ii .le. n) then
         nod = ii         
         if (permut) nod = iperm(nod)          
         if (mask(nod) .eq. maskval) then
c     
c     start a new level
c
            istart = iend
            iend = iend+1 
            riord(iend) = nod
            mask(nod) = 0
            goto 1
         else 
            goto 2
         endif
      endif
c----------------------------------------------------------------------- 
 3    levels(nlev+1) = iend+1 
      do j=1, iend
         mask(riord(j)) = maskval 
      enddo
c----------------------------------------------------------------------- 
      return
      end
c-----------------------------------------------------------------------
      subroutine add_lvst(istart,iend,nlev,riord,ja,ia,mask,maskval) 
      integer nlev, nod, riord(*), ja(*), ia(*), mask(*) 
c-------------------------------------------------------------
c     adds one level set to the previous sets.. 
c     span all nodes of previous mask
c-------------------------------------------------------------
      nod = iend
      do 25 ir = istart+1,iend 
         i = riord(ir)		
         do 24 k=ia(i),ia(i+1)-1
            j = ja(k)
            if (mask(j) .eq. maskval) then
               nod = nod+1 
               mask(j) = 0
               riord(nod) = j
            endif 
 24      continue
 25   continue
      istart = iend 
      iend   = nod 
      return
      end 
c----------------------------------------------------------------------- 
      subroutine stripes (nlev,riord,levels,ip,map,mapptr,ndom)
      implicit none
      integer nlev,riord(*),levels(nlev+1),ip,map(*),
     *    mapptr(*), ndom
c-----------------------------------------------------------------------
c    this is a post processor to BFS. stripes uses the output of BFS to 
c    find a decomposition of the adjacency graph by stripes. It fills 
c    the stripes level by level until a number of nodes .gt. ip is 
c    is reached. 
c---------------------------parameters-----------------------------------
c on entry: 
c --------
c nlev   = number of levels as found by BFS 
c riord  = reverse permutation array produced by BFS -- 
c levels = pointer array for the level structure as computed by BFS. If 
c          lev is a level number, and k1=levels(lev),k2=levels(lev+1)-1, 
c          then all the nodes of level number lev are:
c                      riord(k1),riord(k1+1),...,riord(k2) 
c  ip    = number of desired partitions (subdomains) of about equal size.
c 
c on return
c ---------
c ndom     = number of subgraphs (subdomains) found 
c map      = node per processor list. The nodes are listed contiguously
c            from proc 1 to nproc = mpx*mpy. 
c mapptr   = pointer array for array map. list for proc. i starts at 
c            mapptr(i) and ends at mapptr(i+1)-1 in array map.
c-----------------------------------------------------------------------
c local variables. 
c
      integer ib,ktr,ilev,k,nsiz,psiz 
      ndom = 1 
      ib = 1
c to add: if (ip .le. 1) then ...
      nsiz = levels(nlev+1) - levels(1) 
      psiz = (nsiz-ib)/max(1,(ip - ndom + 1)) + 1 
      mapptr(ndom) = ib 
      ktr = 0 
      do 10 ilev = 1, nlev
c
c     add all nodes of this level to domain
c     
         do 3 k=levels(ilev), levels(ilev+1)-1
            map(ib) = riord(k)
            ib = ib+1
            ktr = ktr + 1 
            if (ktr .ge. psiz  .or. k .ge. nsiz) then 
               ndom = ndom + 1
               mapptr(ndom) = ib 
               psiz = (nsiz-ib)/max(1,(ip - ndom + 1)) + 1 
               ktr = 0
            endif
c
 3       continue
 10   continue
      ndom = ndom-1
      return 
      end
c-----------------------------------------------------------------------
      subroutine stripes0 (ip,nlev,il,ndom,iptr)
      integer ip, nlev, il(*), ndom, iptr(*)
c-----------------------------------------------------------------------
c     This routine is a simple level-set partitioner. It scans
c     the level-sets as produced by BFS from one to nlev.
c     each time the number of nodes in the accumulated set of
c     levels traversed exceeds the parameter ip, this set defines 
c     a new subgraph. 
c-------------------------parameter-list---------------------------------
c on entry:
c --------
c ip     = desired number of nodes per subgraph.
c nlev   = number of levels found  as output by BFS
c il     = integer array containing the pointer array for
c          the level data structure as output by BFS. 
c          thus il(lev+1) - il(lev) = the number of 
c          nodes that constitute the level numbe lev.
c on return
c ---------
c ndom   = number of sungraphs found
c iptr   = pointer array for the sugraph data structure. 
c          thus, iptr(idom) points to the first level that 
c          consistutes the subgraph number idom, in the 
c          level data structure. 
c-----------------------------------------------------------------------
      ktr = 0
      iband = 1 
      iptr(iband) = 1 
c-----------------------------------------------------------------------

      do 10 ilev = 1, nlev
         ktr = ktr + il(ilev+1) - il(ilev)
         if (ktr .gt. ip) then
            iband = iband+1 
            iptr(iband) = ilev+1
            ktr = 0
         endif
c
 10   continue
c-----------returning --------------------
      iptr(iband) = nlev + 1 
      ndom = iband-1
      return
c-----------------------------------------------------------------------
c-----end-of-stripes0--------------------------------------------------- 
      end
c----------------------------------------------------------------------- 
      integer function maskdeg (ja,ia,nod,mask,maskval) 
      implicit none 
      integer ja(*),ia(*),nod,mask(*),maskval
c-----------------------------------------------------------------------
      integer deg, k 
      deg = 0 
      do k =ia(nod),ia(nod+1)-1
         if (mask(ja(k)) .eq. maskval) deg = deg+1 
      enddo
      maskdeg = deg 
      return
      end 
c-----------------------------------------------------------------------
      subroutine perphn(n,ja,ia,init,mask,maskval,nlev,riord,levels) 
      implicit none
      integer n,ja(*),ia(*),init,mask(*),maskval,
     *     nlev,riord(*),levels(*)
c-----------------------------------------------------------------------
c     finds a peripheral node and does a BFS search from it. 
c-----------------------------------------------------------------------
c     see routine  dblstr for description of parameters
c input:
c-------
c ja, ia  = list pointer array for the adjacency graph
c mask    = array used for masking nodes -- see maskval
c maskval = value to be checked against for determing whether or
c           not a node is masked. If mask(k) .ne. maskval then
c           node k is not considered.
c init    = init node in the pseudo-peripheral node algorithm.
c
c output:
c-------
c init    = actual pseudo-peripherial node found.
c nlev    = number of levels in the final BFS traversal.
c riord   =
c levels  =
c----------------------------------------------------------------------- 
      integer j,nlevp,deg,nfirst,mindeg,nod,maskdeg
      integer iperm(1) 
      nlevp = 0 
 1    continue
      riord(1) = init
      nfirst = 1 
      iperm(1) = 0
c
      call BFS(n,ja,ia,nfirst,iperm,mask,maskval,riord,levels,nlev)
      if (nlev .gt. nlevp) then 
         mindeg = n+1 
         do j=levels(nlev),levels(nlev+1)-1
            nod = riord(j) 
            deg = maskdeg(ja,ia,nod,mask,maskval)
            if (deg .lt. mindeg) then
               init = nod
               mindeg = deg
            endif 
         enddo
         nlevp = nlev 
         goto 1 
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine mapper4 (n,ja,ia,ndom,nodes,levst,marker,link)
      implicit none 
      integer n,ndom,ja(*),ia(*),marker(n),levst(2*ndom), 
     *     nodes(*),link(*) 
c-----------------------------------------------------------------------
c     finds domains given ndom centers -- by doing a level set expansion 
c-----------------------------------------------------------------------
c     on entry:
c     ---------
c     n      = dimension of matrix 
c     ja, ia = adajacency list of matrix (CSR format without values) -- 
c     ndom   = number of subdomains (nr output by coarsen)  
c     nodes  = array of size at least n. On input the first ndom entries
c              of nodes should contain the labels of the centers of the 
c              ndom domains from which to do the expansions. 
c     
c     on return 
c     --------- 
c     link  = linked list array for the ndom domains. 
c     nodes = contains the list of nodes of the domain corresponding to
c             link. (nodes(i) and link(i) are related to the same node). 
c    
c     levst = levst(j) points to beginning of subdomain j in link.
c
c     work arrays:
c     ----------- 
c     levst : work array of length 2*ndom -- contains the beginning and
c     end of  current level in link. 
c     beginning of last level in link for each processor.
c     also ends in levst(ndom+i) 
c     marker : work array of length n. 
c
c     Notes on implementation: 
c     ----------------------- 
c     for j .le. ndom link(j) is <0  and indicates the end of the
c     linked list. The most recent element added to the linked
c     list is added at the end of the list (traversal=backward) 
c     For  j .le. ndom, the value of -link(j) is the size of 
c     subdomain j. 
c
c-----------------------------------------------------------------------
c     local variables 
      integer mindom,j,lkend,nod,nodprev,idom,next,i,kk,ii,ilast,nstuck,
     *     isiz, nsize 
c     
c     initilaize nodes and link arrays
c
      do 10 j=1, n
         marker(j) = 0
 10   continue
c     
      do 11 j=1, ndom
         link(j) = -1 
         marker(nodes(j)) = j
         levst(j) = j
         levst(ndom+j) = j 
 11   continue
c
c     ii = next untouched node for restarting new connected component. 
c 
      ii = 0
c     
      lkend = ndom
      nod   = ndom  
      nstuck = 0 
c-----------------------------------------------------------------------
 100  continue 
      idom = mindom(n,ndom,levst,link) 
c-----------------------------------------------------------------------
c     begin level-set loop 
c-----------------------------------------------------------------------
 3    nodprev = nod 
      ilast = levst(ndom+idom) 
      levst(ndom+idom) = lkend       
      next = levst(idom) 
c     
c     linked list traversal loop
c 
      isiz = 0 
      nsize = link(idom) 
 1    i = nodes(next) 
      isiz = isiz + 1 
c     
c     adjacency list traversal loop 
c     
      do 2 kk=ia(i), ia(i+1)-1
         j = ja(kk) 
         if (marker(j) .eq. 0) then 
            call add_lk(j,nod,idom,ndom,lkend,levst,link,nodes,marker) 
         endif
 2    continue
c     
c     if last element of the previous level not reached continue
c     
      if (next .gt. ilast) then
         next = link(next) 
         if (next .gt. 0) goto 1
      endif
c-----------------------------------------------------------------------
c     end level-set traversal --  
c-----------------------------------------------------------------------
      if (nodprev .eq. nod) then
c     
c     link(idom) >0 indicates that set is stuck  --  
c
         link(idom) = -link(idom) 
         nstuck = nstuck+1
      endif
c     
      if (nstuck .lt. ndom) goto 100 
c
c     reset sizes -- 
c 
      do j=1, ndom
         if (link(j) .gt. 0) link(j) = -link(j)
      enddo
c
      if (nod .eq. n) return 
c
c     stuck. add first non assigned point to smallest domain
c     
 20   ii = ii+1
      if (ii .le. n) then
         if (marker(ii) .eq. 0) then 
            idom = 0 
            isiz = n+1 
            do 30 kk=ia(ii), ia(ii+1)-1
               i = marker(ja(kk)) 
               if (i .ne. 0) then 
                  nsize = abs(link(i)) 
                  if (nsize .lt. isiz) then 
                     isiz = nsize
                     idom = i 
                  endif
               endif
 30            continue
c
c     if no neighboring domain select smallest one 
c
               if (idom .eq. 0) idom = mindom(n,ndom,levst,link) 
c     
c     add ii to sudomain idom at end of linked list  
c     
            call add_lk(ii,nod,idom,ndom,lkend,levst,link,nodes,marker) 
            goto 3 
         else
            goto 20 
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine get_domns2(ndom,nodes,link,levst,riord,iptr)
      implicit none 
      integer ndom,nodes(*),link(*),levst(*),riord(*),iptr(*) 
c-----------------------------------------------------------------------
c     constructs the subdomains from its linked list data structure
c-----------------------------------------------------------------------
c     input:
c     ndom  = number of subdomains
c     nodes = sequence of nodes are produced by mapper4. 
c     link  = link list array as produced by mapper4.
c     on return:
c----------
c     riord = contains the nodes in each subdomain in succession.
c     iptr  = pointer in riord for beginnning of each subdomain.
c     Thus subdomain number i consists of nodes 
c     riord(k1),riord(k1)+1,...,riord(k2) 
c     where k1 = iptr(i), k2= iptr(i+1)-1
c     
c-----------------------------------------------------------------------
c     local variables 
      integer nod, j, next, ii 
      nod = 1
      iptr(1) = nod 
      do 21 j=1, ndom 
         next = levst(j)
 22      ii = nodes(next)
         riord(nod) = ii 
         nod = nod+1 
         next = link(next) 
         if (next .gt.  0) goto 22
         iptr(j+1) = nod 
 21   continue
c
      return
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      function mindom(n, ndom, levst, link) 
      implicit none 
      integer mindom, n, ndom, levst(2*ndom),link(n) 
c-----------------------------------------------------------------------
c     returns  the domain with smallest size
c----------------------------------------------------------------------- 
c      locals
c
      integer i, nsize, isiz 
c
      isiz = n+1 
      do 10 i=1, ndom
         nsize = - link(i) 
         if (nsize .lt. 0) goto 10 
         if (nsize .lt. isiz) then 
            isiz = nsize
            mindom = i
         endif
 10   continue
      return
      end 
c-----------------------------------------------------------------------
      subroutine add_lk(new,nod,idom,ndom,lkend,levst,link,nodes,marker) 
      implicit none
      integer new,nod,idom,ndom,lkend,levst(*),link(*),nodes(*),
     *     marker(*) 
c----------------------------------------------------------------------- 
c     adds from head -- 
c
c     adds one entry (new) to linked list and ipdates everything.
c     new  = node to be added
c     nod  = current number of marked nodes
c     idom = domain to which new is to be added
c     ndom = total number of domains
c     lkend= location of end of structure (link and nodes)
c     levst= pointer array for link, nodes
c     link = link array 
c     nodes= nodes array -- 
c     marker = marker array == if marker(k) =0 then node k is not
c              assigned yet. 
c----------------------------------------------------------------------- 
c      locals
c     
      integer ktop  
      lkend = lkend + 1
      nodes(lkend) = new
      nod = nod+1 
      marker(new) = idom 
      ktop = levst(idom) 
      link(lkend) = ktop 
      link(idom) = link(idom)-1 
      levst(idom) = lkend 
      return
c-----------------------------------------------------------------------
c-------end-of-add_lk--------------------------------------------------- 
      end 
c----------------------------------------------------------------------- 
      subroutine find_ctr(n,nsiz,ja,ia,init,mask,maskval,riord,
     *     levels,center,iwk) 
      implicit none
      integer n,nsiz,ja(*),ia(*),init,mask(*),maskval,riord(*),
     *     levels(*),center,iwk(*) 
c-----------------------------------------------------------------------
c     finds a center point of a subgraph -- 
c-----------------------------------------------------------------------
c     n, ja, ia = graph
c     nsiz = size of current domain.
c     init = initial node in search
c     mask
c     maskval 
c-----------------------------------------------------------------------
c     local variables 
      integer midlev, nlev,newmask, k, kr, kl, init0, nlev0  
      call perphn(n,ja,ia,init,mask,maskval,nlev,riord,levels)
c-----------------------------------------------------------------------
c     midlevel = level which cuts domain into 2 roughly equal-size 
c     regions 
c
      midlev = 1
      k = 0 
 1    continue
      k = k + levels(midlev+1)-levels(midlev) 
      if (k*2 .lt. nsiz) then
         midlev = midlev+1
         goto 1 
      endif
c-----------------------------------------------------------------------
      newmask = n+maskval
c     
c     assign temporary masks to mid-level elements
c     
      do k=levels(midlev),levels(midlev+1)-1
         mask(riord(k)) = newmask
      enddo
c     
c     find pseudo-periph node for mid-level `line'
c
      kr = 1
      kl = kr + nsiz 
      init0 = riord(levels(midlev))
      call perphn(n,ja,ia,init0,mask,newmask,nlev0,iwk(kr),iwk(kl)) 
c-----------------------------------------------------------------------
c     restore  mask to initial state 
c-----------------------------------------------------------------------
      do k=levels(midlev),levels(midlev+1)-1
         mask(riord(k)) = maskval 
      enddo
c-----------------------------------------------------------------------
c     define center 
c-----------------------------------------------------------------------  
      midlev = 1 + (nlev0-1)/2
      k = iwk(kl+midlev-1)
      center = iwk(k) 
c----------------------------------------------------------------------- 
      return 
      end 
c-----------------------------------------------------------------------
      subroutine rversp (n, riord)
      integer n, riord(n)
c-----------------------------------------------------------------------
c     this routine does an in-place reversing of the permutation array
c     riord --
c-----------------------------------------------------------------------
      integer j, k
      do 26 j=1,n/2
         k = riord(j)
         riord(j) = riord(n-j+1)
         riord(n-j+1) = k
 26   continue
      return
      end
c-----------------------------------------------------------------------
