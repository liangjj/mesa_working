      subroutine cgevlrv(a,lda,n,e,lv,rv,ldv,work,task,info)
c***begin prologue  cgeev
c***date written   800808   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d4a4
c***keywords  complex,eigenvalue,eigenvector,general matrix
c***author  kahaner, d. k., (nbs)
c           moler, c. b., (u. of new mexico)
c           stewart, g. w., (u. of maryland)
c***purpose  to compute the eigenvalues and, optionally, the eigen-
c            vectors of a general complex matrix.
c***description
c
c     licepack.    this version dated 08/08/80.
c     david kahner, cleve moler, g. w. stewart
c       n.b.s.         u.n.m.      n.b.s./u.md.
c
c     abstract
c      cgeev computes the eigenvalues and, optionally,
c      the eigenvectors of a general complex matrix.
c
c     call sequence parameters-
c       (the values of parameters marked with * (star) will be changed
c         by cgeev.)
c
c        a*      complex(lda,n)
c                complex nonsymmetric input matrix.
c
c        lda     integer
c                set by the user to
c                the leading dimension of the complex array a.
c
c        n       integer
c                set by the user to
c                the order of the matrices a and v, and
c                the number of elements in e.
c
c        e*      complex(n)
c                on return from cgeev e contains the eigenvalues of a.
c                see also task below.
c
c        lv*     complex(ldv,n)
c                if the user has set task 
c                = eigenvalues only  lv is not referenced.
c                = left eigenvectors  the n left eigenvectors of a 
c                  are stored in the first n columns of lv.
c                (if the input matrix a is nearly degenerate, lv
c                 will be badly conditioned, i.e. have nearly
c                 dependent columns.)
c
c        rv*     complex(ldv,n)
c                if the user has set task
c                = eigenvalues only  rv is not referenced.
c                = right eigenvectors  the n right eigenvectors of a are 
c                  stored in the first n columns of rv.
c                (if the input matrix a is nearly degenerate, v
c                 will be badly conditioned, i.e. have nearly
c                 dependent columns.)
c
c        ldv     integer
c                set by the user to
c                the leading dimension of the array lv and rv.
c                n must be .le. ldv.
c
c        work*   real(3n)
c                temporary storage vector.  contents changed by cgeev.
c
c        task    character
c                set by the user to
c                eigenvalues only.
c                right eigenvectors            eigenvalues and right vectors
c                                              to be calculated.
c                left eigenvectors             eigenvalues and left vectors
c                                              to be calculated.
c                left and right eigenvectors   eigenvalues and both left and
c                                              right vectors calculated.
c
c        info*   integer
c                on return from cgeev the value of info is
c                = 0  normal return, calculation successful.
c                = k  if the eigenvalue iteration fails to converge,
c                     eigenvalues k+1 through n are correct, but
c                     no eigenvectors were computed even if they were
c                     requested (job nonzero).
c
c      error messages
c           no. 1  recoverable  n is greater than lda
c           no. 2  recoverable  n is less than one.
c           no. 3  recoverable  job is nonzero and n is greater than ldv
c           no. 4  warning      lda > ldv,  elements of a other than the
c          n by n input elements have been changed
c           no. 5  warning      lda < ldv,  elements of v other than the
c          n by n output elements have been changed
c
c
c     subroutines used
c
c     eispack-  cbabk2, cbal, comqr, comqr2, corth
c     blas-  scopy
c     mesalib- lnkerr
c***references  (none)
c***routines called  cbabk2,cbal,comqr,comqr2,corth,scopy,lnkerr
c***end prologue  cgeev
      implicit real *8 (a-h,o-z)
      real *8 a(*), lv(*), rv(*), e(*), work(*), dum
      character*(*) task
      common/io/inp,iout
c***first executable statement  cgeev
      if(n .gt. lda) call lnkerr( 'cgeev-n .gt. lda')
      if(n .gt. lda) return
      if(n .lt. 1) call lnkerr( 'cgeev-n .lt. 1')
      if(n .lt. 1) return
      if(n .eq. 1 .and. task .eq. 'eigenvalues only') go to 35
      mdim = 2 * lda
      if(task .eq. 'eigenvalues only') go to 5
      if(n .gt. ldv) call lnkerr( 'cgeev-job .ne. 0, and n .gt. ldv')
      if(n .gt. ldv) return
      if(n .eq. 1) go to 35
c
c       rearrange a if necessary when lda.gt.ldv and task
c       .ne.'eigenvalues only'
c
      mdim = min(mdim,2 * ldv)
      if(lda.lt.ldv) call lnkerr( 'cgeev-lda.lt.ldv,  elements of v othe
     1r than the n by n output elements have been changed.')
      if(lda.le.ldv) go to 5
      call lnkerr(  'cgeev-lda.gt.ldv, elements of a other than the n by
     1 n input elements have been changed.')
      l = n - 1
      do 4 j=1,l
         i = 2 * n
         m = 1+j*2*ldv
         k = 1+j*2*lda
         call scopy(i,a(k),1,a(m),1)
    4 continue
    5 continue
c
c     separate real and imaginary parts
c
      do 6 j = 1, n
       k = (j-1) * mdim +1
       l = k + n
       call scopy(n,a(k+1),2,work(1),1)
       call scopy(n,a(k),2,a(k),1)
       call scopy(n,work(1),1,a(l),1)
    6 continue
c
c     scale and orthogonal reduction to hessenberg.
c
      call cbal(mdim,n,a(1),a(n+1),ilo,ihi,work(1))
      call corth(mdim,n,ilo,ihi,a(1),a(n+1),work(n+1),work(2*n+1))
c
      if(task .ne. 'eigenvalues only') go to 10
c
c     eigenvalues only
c
      call comqr(mdim,n,ilo,ihi,a(1),a(n+1),e(1),e(n+1),info)
      go to 30
c
c     eigenvalues and eigenvectors.
c
   10 continue
c
c     at this point we can compute left, right or left and right
c     eigenvectors.
c
c     in the routine below we calculate the eigenvalues and eigenvectors
c     of the upper hessenberg matrix.  they are stored in either vl or vr
c     depending on task.
c
      if(task.eq.'right eigenvectors') then
         call comqr2(mdim,n,ilo,ihi,work(n+1),work(2*n+1),a(1),a(n+1),
     1               e(1),e(n+1),rv(1),rv(n+1),info)
         if (info .ne. 0) go to 30
         call cbabk3(mdim,n,ilo,ihi,work(1),n,rv(1),rv(n+1),dum,dum,
     1               work(n+1),task)
c
c     convert eigenvectors to complex storage.
c
      do 20 j = 1,n
       k = (j-1) * mdim + 1
       i = (j-1) * 2 * ldv + 1
       l = k + n
       call scopy(n,rv(k),1,work(1),1)
       call scopy(n,rv(l),1,rv(i+1),2)
       call scopy(n,work(1),1,rv(i),2)
   20 continue
      elseif(task.eq.'left eigenvectors') then
         call comqr2(mdim,n,ilo,ihi,work(n+1),work(2*n+1),a(1),a(n+1),
     1               e(1),e(n+1),lv(1),lv(n+1),info)
         if (info .ne. 0) go to 30
         call cbabk3(mdim,n,ilo,ihi,work(1),n,dum,dum,lv(1),lv(n+1),
     1               work(n+1),task)
c
c     convert eigenvectors to complex storage.
c
      do 200 j = 1,n
       k = (j-1) * mdim + 1
       i = (j-1) * 2 * ldv + 1
       l = k + n
       call scopy(n,lv(k),1,work(1),1)
       call scopy(n,lv(l),1,lv(i+1),2)
       call scopy(n,work(1),1,lv(i),2)
 200  continue
      elseif(task.eq.'left and right eigenvectors') then
         call comqr2(mdim,n,ilo,ihi,work(n+1),work(2*n+1),a(1),a(n+1),
     1               e(1),e(n+1),rv(1),rv(n+1),info)
         if (info .ne. 0) go to 30
         call cbabk3(mdim,n,ilo,ihi,work(1),n,rv(1),rv(n+1),
     1               lv(1),lv(n+1),work(n+1),task)
c
c     convert eigenvectors to complex storage.
c
      do 300 j = 1,n
       k = (j-1) * mdim + 1
       i = (j-1) * 2 * ldv + 1
       l = k + n
       call scopy(n,rv(k),1,work(1),1)
       call scopy(n,rv(l),1,rv(i+1),2)
       call scopy(n,work(1),1,rv(i),2)
       call scopy(n,lv(k),1,work(1),1)
       call scopy(n,lv(l),1,lv(i+1),2)
       call scopy(n,work(1),1,lv(i),2)
 300  continue
      endif
c
c     convert eigenvalues to complex storage.
c
   30 call scopy(n,e(1),1,work(1),1)
      call scopy(n,e(n+1),1,e(2),2)
      call scopy(n,work(1),1,e(1),2)
      return
c
c     take care of n=1 case
c
   35 e(1) = a(1)
      e(2) = a(2)
      info = 0
      if(task .eq. 'eigenvalues only') then
         return
      endif
      if(task.eq.'right eigenvectors') then
         rv(1) = 1.d0
         rv(2) = 0.d0
      elseif(task.eq.'left eigenvectors') then
         lv(1) = 1.d0
         lv(2) = 0.d0
      elseif(task.eq.'left and right eigenvectors') then
         rv(1) = 1.d0
         rv(2) = 0.d0
         lv(1) = 1.d0
         lv(2) = 0.d0
      endif
      return
      end
