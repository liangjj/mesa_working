*deck sspev
      subroutine sspev (a, n, e, v, ldv, work, job, info)
c***begin prologue  sspev
c***purpose  compute the eigenvalues and, optionally, the eigenvectors
c            of a real symmetric matrix stored in packed form.
c***library   slatec (eispack)
c***category  d4a1
c***type      single precision (sspev-s)
c***keywords  eigenvalues, eigenvectors, eispack, packed, symmetric
c***author  kahaner, d. k., (nbs)
c           moler, c. b., (u. of new mexico)
c           stewart, g. w., (u. of maryland)
c***description
c
c     abstract
c      sspev computes the eigenvalues and, optionally, the eigenvectors
c      of a real symmetric matrix stored in packed form.
c
c     call sequence parameters-
c       (the values of parameters marked with * (star) will be  changed
c         by sspev.)
c
c        a*      real(n*(n+1)/2)
c                real symmetric packed input matrix.  contains upper
c                triangle and diagonal of a, by column (elements
c                11, 12, 22, 13, 23, 33, ...).
c
c        n       integer
c                set by the user to
c                the order of the matrix a.
c
c        e*      real(n)
c                on return from sspev, e contains the eigenvalues of a.
c                see also info below.
c
c        v*      real(ldv,n)
c                on return from sspev, if the user has set job
c                = 0        v is not referenced.
c                = nonzero  the n eigenvectors of a are stored in the
c                first n columns of v.  see also info below.
c
c        ldv     integer
c                set by the user to
c                the leading dimension of the array v if job is also
c                set nonzero.  in that case, n must be .le. ldv.
c                if job is set to zero, ldv is not referenced.
c
c        work*   real(2n)
c                temporary storage vector.  contents changed by sspev.
c
c        job     integer
c                set by the user to
c                = 0        eigenvalues only to be calculated by sspev.
c                           neither v nor ldv are referenced.
c                = nonzero  eigenvalues and vectors to be calculated.
c                           in this case, a & v must be distinct arrays.
c                           also, if lda .gt. ldv, sspev changes all the
c                           elements of a thru column n.  if lda < ldv,
c                           sspev changes all the elements of v through
c                           column n.  if lda=ldv, only a(i,j) and v(i,
c                           j) for i,j = 1,...,n are changed by sspev.
c
c       info*   integer
c               on return from sspev, the value of info is
c               = 0 for normal return.
c               = k if the eigenvalue iteration fails to converge.
c                   eigenvalues and vectors 1 through k-1 are correct.
c
c
c     error messages-
c          no. 1   recoverable  n is greater than ldv and job is nonzero
c          no. 2   recoverable  n is less than one
c
c***references  (none)
c***routines called  imtql2, tqlrat, trbak3, tred3, xermsg
c***revision history  (yymmdd)
c   800808  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  sspev
      integer i,info,j,ldv,m,n
      real a(*),e(*),v(ldv,*),work(*)
c***first executable statement  sspev
       if (n .gt. ldv) call xermsg ('slatec', 'sspev', 'n .gt. ldv.',
     +    1, 1)
       if(n .gt. ldv) return
       if (n .lt. 1) call xermsg ('slatec', 'sspev', 'n .lt. 1', 2, 1)
       if(n .lt. 1) return
c
c       check n=1 case
c
      e(1) = a(1)
      info = 0
      if(n .eq. 1) return
c
      if(job.ne.0) go to 20
c
c     eigenvalues only
c
      call tred3(n,1,a,e,work(1),work(n+1))
      call tqlrat(n,e,work(n+1),info)
      return
c
c     eigenvalues and eigenvectors
c
   20 call tred3(n,1,a,e,work(1),work(1))
      do 30 i = 1, n
        do 25 j = 1, n
   25     v(i,j) = 0.
   30   v(i,i) = 1.
      call imtql2(ldv,n,e,work,v,info)
      m = n
      if(info .ne. 0) m = info - 1
      call trbak3(ldv,n,1,a,m,v)
      return
      end
