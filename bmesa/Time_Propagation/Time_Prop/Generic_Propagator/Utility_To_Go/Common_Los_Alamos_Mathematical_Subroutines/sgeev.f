*deck sgeev
      subroutine sgeev (a, lda, n, e, v, ldv, work, job, info)
c***begin prologue  sgeev
c***purpose  compute the eigenvalues and, optionally, the eigenvectors
c            of a real general matrix.
c***library   slatec
c***category  d4a2
c***type      single precision (sgeev-s, cgeev-c)
c***keywords  eigenvalues, eigenvectors, general matrix
c***author  kahaner, d. k., (nbs)
c           moler, c. b., (u. of new mexico)
c           stewart, g. w., (u. of maryland)
c***description
c
c     abstract
c      sgeev computes the eigenvalues and, optionally,
c      the eigenvectors of a general real matrix.
c
c     call sequence parameters-
c       (the values of parameters marked with * (star) will be changed
c         by sgeev.)
c
c        a*      real(lda,n)
c                real nonsymmetric input matrix.
c
c        lda     integer
c                set by the user to
c                the leading dimension of the real array a.
c
c        n       integer
c                set by the user to
c                the order of the matrices a and v, and
c                the number of elements in e.
c
c        e*      complex(n)
c                on return from sgeev, e contains the eigenvalues of a.
c                see also info below.
c
c        v*      complex(ldv,n)
c                on return from sgeev, if the user has set job
c                = 0        v is not referenced.
c                = nonzero  the n eigenvectors of a are stored in the
c                first n columns of v.  see also info below.
c                (note that if the input matrix a is nearly degenerate,
c                 v may be badly conditioned, i.e., may have nearly
c                 dependent columns.)
c
c        ldv     integer
c                set by the user to
c                the leading dimension of the array v if job is also
c                set nonzero.  in that case, n must be .le. ldv.
c                if job is set to zero, ldv is not referenced.
c
c        work*   real(2n)
c                temporary storage vector.  contents changed by sgeev.
c
c        job     integer
c                set by the user to
c                = 0        eigenvalues only to be calculated by sgeev.
c                           neither v nor ldv is referenced.
c                = nonzero  eigenvalues and vectors to be calculated.
c                           in this case, a & v must be distinct arrays.
c                           also, if lda .gt. ldv, sgeev changes all the
c                           elements of a thru column n.  if lda < ldv,
c                           sgeev changes all the elements of v through
c                           column n. if lda = ldv, only a(i,j) and v(i,
c                           j) for i,j = 1,...,n are changed by sgeev.
c
c        info*   integer
c                on return from sgeev the value of info is
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
c           no. 4  warning      lda > ldv, elements of a other than the
c                               n by n input elements have been changed.
c           no. 5  warning      lda < ldv, elements of v other than the
c                               n x n output elements have been changed.
c
c***references  (none)
c***routines called  balanc, balbak, hqr, hqr2, orthes, ortran, scopy,
c                    scopym, xermsg
c***revision history  (yymmdd)
c   800808  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  sgeev
      integer i,ihi,ilo,info,j,jb,job,k,km,kp,l,lda,ldv,
     1        mdim,n
      real*8 a(*),e(*),work(*),v(*)
c***first executable statement  sgeev
      if (n .gt. lda) call lnkerr ('sgeev n .gt. lda.')
      if (n .gt. lda) return
      if (n .lt. 1) call lnkerr ('sgeev-n .lt. 1')
      if(n .lt. 1) return
      if(n .eq. 1 .and. job .eq. 0) go to 35
      mdim = lda
      if(job .eq. 0) go to 5
      if (n .gt. ldv) call lnkerr ('sgeev-job .ne. 0, and n .gt. ldv')
      if(n .gt. ldv) return
      if(n .eq. 1) go to 35
c
c       rearrange a if necessary when lda.gt.ldv and job .ne.0
c
      mdim = min(lda,ldv)
      if(lda.lt.ldv) call lnkerr( 'sgeev-lda.lt.ldv,  elements of v othe
     1r than the n by n output elements have been changed.')
      if(lda.le.ldv) go to 5
      call lnkerr(  'sgeev-lda.gt.ldv, elements of a other than the n by
     1 n input elements have been changed.')
      l = n - 1
      do 4 j=1,l
         m = 1+j*ldv
         k = 1+j*lda
         call scopy(n,a(k),1,a(m),1)
    4 continue
    5 continue
c
c     scale and orthogonal reduction to hessenberg.
c
      call balanc(mdim,n,a,ilo,ihi,work(1))
      call orthes(mdim,n,ilo,ihi,a,work(n+1))
      if(job .ne. 0) go to 10
c
c     eigenvalues only
c
      call hqr(lda,n,ilo,ihi,a,e(1),e(n+1),info)
      go to 30
c
c     eigenvalues and eigenvectors.
c
   10 call ortran(mdim,n,ilo,ihi,a,work(n+1),v)
      call hqr2(mdim,n,ilo,ihi,a,e(1),e(n+1),v,info)
      if (info .ne. 0) go to 30
      call balbak(mdim,n,ilo,ihi,work(1),n,v)
c
c     convert eigenvectors to complex storage.
c
      do 20 jb = 1,n
         j=n+1-jb
         i=n+j
         k=(j-1)*mdim+1
         kp=k+mdim
         km=k-mdim
         if(e(i).ge.0.0d0) call scopy(n,v(k),1,work(1),2)
         if(e(i).lt.0.0d0) call scopy(n,v(km),1,work(1),2)
         if(e(i).eq.0.0d0) call scopy(n,0.0d0,0,work(2),2)
         if(e(i).gt.0.0d0) call scopy(n,v(kp),1,work(2),2)
         if(e(i).lt.0.0d0) call scopym(n,v(k),1,work(2),2)
         l=2*(j-1)*ldv+1
         call scopy(2*n,work(1),1,v(l),1)
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
      e(2) = 0.d0
      info = 0
      if(job .eq. 0) return
      v(1) = a(1)
      v(2) = 0.d0
      return
      end
