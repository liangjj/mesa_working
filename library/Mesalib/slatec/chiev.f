*deck chiev
      subroutine chiev (a, lda, n, e, v, ldv, work, job, info)
c***begin prologue  chiev
c***purpose  compute the eigenvalues and, optionally, the eigenvectors
c            of a complex hermitian matrix.
c***library   slatec
c***category  d4a3
c***type      complex (ssiev-s, chiev-c)
c***keywords  complex hermitian, eigenvalues, eigenvectors, matrix,
c             symmetric
c***author  kahaner, d. k., (nbs)
c           moler, c. b., (u. of new mexico)
c           stewart, g. w., (u. of maryland)
c***description
c
c     david kahaner, cleve moler, g. w. stewart,
c       n.b.s.         u.n.m.      n.b.s./u.md.
c
c     abstract
c      chiev computes the eigenvalues and, optionally,
c      the eigenvectors of a complex hermitian matrix.
c
c     call sequence parameters-
c       (the values of parameters marked with * (star) will be changed
c         by chiev.)
c
c        a*      complex(lda,n)
c                complex hermitian input matrix.
c                only the upper triangle of a need be
c                filled in.  elements on diagonal must be real.
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
c        e*      real(n)
c                on return from chiev e contains the eigenvalues of a.
c                see also info below.
c
c        v*      complex(ldv,n)
c                on return from chiev if the user has set job
c                = 0        v is not referenced.
c                = nonzero  the n eigenvectors of a are stored in the
c                first n columns of v.  see also info below.
c
c        ldv     integer
c                set by the user to
c                the leading dimension of the array v if job is also
c                set nonzero.  in that case n must be .le. ldv.
c                if job is set to zero ldv is not referenced.
c
c        work*   real(4n)
c                temporary storage vector.  contents changed by chiev.
c
c        job     integer
c                set by the user to
c                = 0        eigenvalues only to be calculated by chiev.
c                           neither v nor ldv are referenced.
c                = nonzero  eigenvalues and vectors to be calculated.
c                           in this case a and v must be distinct arrays
c                           also if lda .gt. ldv chiev changes all the
c                           elements of a thru column n.  if lda < ldv
c                           chiev changes all the elements of v through
c                           column n.  if lda = ldv only a(i,j) and v(i,
c                           j) for i,j = 1,...,n are changed by chiev.
c
c        info*   integer
c                on return from chiev the value of info is
c                = 0  normal return, calculation successful.
c                = k  if the eigenvalue iteration fails to converge,
c                     eigenvalues (and eigenvectors if requested)
c                     1 through k-1 are correct.
c
c      error messages
c           no. 1  recoverable  n is greater than lda
c           no. 2  recoverable  n is less than one.
c           no. 3  recoverable  job is nonzero and n is greater than ldv
c           no. 4  warning      lda > ldv,  elements of a other than the
c                               n by n input elements have been changed
c           no. 5  warning      lda < ldv,  elements of v other than the
c                               n by n output elements have been changed
c           no. 6  recoverable  nonreal element on diagonal of a.
c
c***references  (none)
c***routines called  htribk, htridi, imtql2, scopy, scopym, tqlrat,
c                    xermsg
c***revision history  (yymmdd)
c   800808  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  chiev
      integer i,info,j,job,k,l,lda,ldv,m,mdim,n
      real a(*),e(*),work(*),v(*)
c***first executable statement  chiev
      if (n .gt. lda) call xermsg ('slatec', 'chiev', 'n .gt. lda.', 1,
     +   1)
      if(n .gt. lda) return
      if (n .lt. 1) call xermsg ('slatec', 'chiev', 'n .lt. 1', 2, 1)
      if(n .lt. 1) return
      if(n .eq. 1 .and. job .eq. 0) go to 35
      mdim = 2 * lda
      if(job .eq. 0) go to 5
      if (n .gt. ldv) call xermsg ('slatec', 'chiev',
     +   'job .ne. 0 and n .gt. ldv.', 3, 1)
      if(n .gt. ldv) return
      if(n .eq. 1) go to 35
c
c       rearrange a if necessary when lda.gt.ldv and job .ne.0
c
      mdim = min(mdim,2 * ldv)
      if (lda .lt. ldv) call xermsg ('slatec', 'chiev',
     +   'lda.lt.ldv,  elements of v other than the n by n output ' //
     +   'elements have been changed.', 5, 0)
      if(lda.le.ldv) go to 5
      call xermsg ('slatec', 'chiev',
     +   'lda.gt.ldv, elements of a other than the n by n input ' //
     +   'elements have been changed.', 4, 0)
      l = n - 1
      do 4 j=1,l
         m = 1+j*2*ldv
         k = 1+j*2*lda
         call scopy(2*n,a(k),1,a(m),1)
    4 continue
    5 continue
c
c     fill in lower triangle of a, column by column.
c
      do 6 j = 1,n
       k = (j-1)*(mdim+2)+1
       if (a(k+1) .ne. 0.0) call xermsg ('slatec', 'chiev',
     +    'nonreal element on diagonal of a', 6, 1)
      if(a(k+1) .ne.0.0) return
       call scopy(n-j+1,a(k),mdim,a(k),2)
       call scopym(n-j+1,a(k+1),mdim,a(k+1),2)
    6 continue
c
c     separate real and imaginary parts
c
      do 10 j = 1, n
       k = (j-1) * mdim +1
       l = k + n
       call scopy(n,a(k+1),2,work(1),1)
       call scopy(n,a(k),2,a(k),1)
       call scopy(n,work(1),1,a(l),1)
   10 continue
c
c    reduce a to tridiagonal matrix.
c
      call htridi(mdim,n,a(1),a(n+1),e,work(1),work(n+1),
     1            work(2*n+1))
      if(job .ne. 0) goto 15
c
c     eigenvalues only.
c
      call tqlrat(n,e,work(n+1),info)
      return
c
c     eigenvalues and eigenvectors.
c
   15 do 17 j = 1,n
       k = (j-1) * mdim + 1
       m = k + n - 1
       do 16 i = k,m
   16   v(i) = 0.
       i = k + j - 1
       v(i) = 1.
   17 continue
      call imtql2(mdim,n,e,work(1),v,info)
      if(info .ne. 0) return
      call htribk(mdim,n,a(1),a(n+1),work(2*n+1),n,v(1),v(n+1))
c
c    convert eigenvectors to complex storage.
c
      do 20 j = 1,n
       k = (j-1) * mdim + 1
       i = (j-1) * 2 * ldv + 1
       l = k + n
       call scopy(n,v(k),1,work(1),1)
       call scopy(n,v(l),1,v(i+1),2)
       call scopy(n,work(1),1,v(i),2)
   20 continue
      return
c
c     take care of n=1 case.
c
   35 if (a(2) .ne. 0.) call xermsg ('slatec', 'chiev',
     +   'nonreal element on diagonal of a', 6, 1)
      if(a(2) .ne. 0.) return
      e(1) = a(1)
      info = 0
      if(job .eq. 0) return
      v(1) = a(1)
      v(2) = 0.
      return
      end
