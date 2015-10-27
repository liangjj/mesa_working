      subroutine cgeev(a,lda,n,e,v,ldv,work,job,info)
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
c                see also info below.
c
c        v*      complex(ldv,n)
c                on return from cgeev if the user has set job
c                = 0        v is not referenced.
c                = nonzero  the n eigenvectors of a are stored in the
c                first n columns of v.  see also info below.
c                (if the input matrix a is nearly degenerate, v
c                 will be badly conditioned, i.e. have nearly
c                 dependent columns.)
c
c        ldv     integer
c                set by the user to
c                the leading dimension of the array v if job is also
c                set nonzero.  in that case n must be .le. ldv.
c                if job is set to zero ldv is not referenced.
c
c        work*   real(3n)
c                temporary storage vector.  contents changed by cgeev.
c
c        job     integer
c                set by the user to
c                = 0        eigenvalues only to be calculated by cgeev.
c                           neither v nor ldv are referenced.
c                = nonzero  eigenvalues and vectors to be calculated.
c                           in this case a & v must be distinct arrays.
c                           also,  if lda > ldv,  cgeev changes all the
c                           elements of a thru column n.  if lda < ldv,
c                           cgeev changes all the elements of v through
c                           column n.  if lda = ldv only a(i,j) and v(i,
c                           j) for i,j = 1,...,n are changed by cgeev.
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
      real *8 a(*), v(*), e(*), work(*)
      common/io/inp,iout
c***first executable statement  cgeev
      if(n .gt. lda) call lnkerr( 'cgeev-n .gt. lda')
      if(n .gt. lda) return
      if(n .lt. 1) call lnkerr( 'cgeev-n .lt. 1')
      if(n .lt. 1) return
      if(n .eq. 1 .and. job .eq. 0) go to 35
      mdim = 2 * lda
      if(job .eq. 0) go to 5
      if(n .gt. ldv) call lnkerr( 'cgeev-job .ne. 0, and n .gt. ldv')
      if(n .gt. ldv) return
      if(n .eq. 1) go to 35
c
c       rearrange a if necessary when lda.gt.ldv and job .ne.0
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
      if(job .ne. 0) go to 10
c
c     eigenvalues only
c
      call comqr(mdim,n,ilo,ihi,a(1),a(n+1),e(1),e(n+1),info)
      go to 30
c
c     eigenvalues and eigenvectors.
c
   10 continue
      call comqr2(mdim,n,ilo,ihi,work(n+1),work(2*n+1),a(1),a(n+1),
     1  e(1),e(n+1),v(1),v(n+1),info)
      if (info .ne. 0) go to 30
      call cbabk2(mdim,n,ilo,ihi,work(1),n,v(1),v(n+1))
c
c     convert eigenvectors to complex storage.
c
      do 20 j = 1,n
       k = (j-1) * mdim + 1
       i = (j-1) * 2 * ldv + 1
       l = k + n
       call scopy(n,v(k),1,work(1),1)
       call scopy(n,v(l),1,v(i+1),2)
       call scopy(n,work(1),1,v(i),2)
   20 continue
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
      if(job .eq. 0) return
c      v(1) = a(1)
c      v(2) = a(2)
      v(1)=1.d0
      v(2)=0.d0
      return
      end
