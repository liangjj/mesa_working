*deck ssiev
      subroutine ssiev (a, lda, n, e, work, job, info)
c***begin prologue  ssiev
c***purpose  compute the eigenvalues and, optionally, the eigenvectors
c            of a real symmetric matrix.
c***library   slatec
c***category  d4a1
c***type      single precision (ssiev-s, chiev-c)
c***keywords  complex hermitian, eigenvalues, eigenvectors, matrix,
c             symmetric
c***author  kahaner, d. k., (nbs)
c           moler, c. b., (u. of new mexico)
c           stewart, g. w., (u. of maryland)
c***description
c
c     abstract
c      ssiev computes the eigenvalues and, optionally, the eigenvectors
c      of a real symmetric matrix.
c
c     call sequence parameters-
c       (the values of parameters marked with * (star) will be  changed
c         by ssiev.)
c
c       a*      real (lda,n)
c               real symmetric input matrix.
c               only the diagonal and upper triangle of a must be input,
c               as ssiev copies the upper triangle to the lower.
c               that is, the user must define a(i,j), i=1,..n, and j=i,.
c               ..,n.
c               on return from ssiev, if the user has set job
c               = 0        the lower triangle of a has been altered.
c               = nonzero  the n eigenvectors of a are stored in its
c               first n columns.  see also info below.
c
c       lda     integer
c               set by the user to
c               the leading dimension of the array a.
c
c       n       integer
c               set by the user to
c               the order of the matrix a and
c               the number of elements in e.
c
c       e*      real (n)
c               on return from ssiev, e contains the n
c               eigenvalues of a.  see also info below.
c
c       work*   real (2*n)
c               temporary storage vector.  contents changed by ssiev.
c
c       job     integer
c               set by user on input
c               = 0         only calculate eigenvalues of a.
c               = nonzero   calculate eigenvalues and eigenvectors of a.
c
c       info*   integer
c               on return from ssiev, the value of info is
c               = 0 for normal return.
c               = k if the eigenvalue iteration fails to converge.
c                   eigenvalues and vectors 1 through k-1 are correct.
c
c
c     error messages-
c          no. 1   recoverable  n is greater than lda
c          no. 2   recoverable  n is less than one
c
c***references  (none)
c***routines called  imtql2, tqlrat, tred1, tred2, xermsg
c***revision history  (yymmdd)
c   800808  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  ssiev
      integer info,job,lda,n
      real a(lda,*),e(*),work(*)
c***first executable statement  ssiev
       if (n .gt. lda) call xermsg ('slatec', 'ssiev', 'n .gt. lda.',
     +       1, 1)
       if(n .gt. lda) return
       if (n .lt. 1) call xermsg ('slatec', 'ssiev', 'n .lt. 1', 2, 1)
       if(n .lt. 1) return
c
c       check n=1 case
c
      e(1) = a(1,1)
      info = 0
      if(n .eq. 1) return
c
c     copy upper triangle to lower
c
      do 10 j=1,n
      do 10 i=1,j
         a(j,i)=a(i,j)
   10 continue
c
      if(job.ne.0) go to 20
c
c     eigenvalues only
c
      call tred1(lda,n,a,e,work(1),work(n+1))
      call tqlrat(n,e,work(n+1),info)
      return
c
c     eigenvalues and eigenvectors
c
   20 call tred2(lda,n,a,e,work,a)
      call imtql2(lda,n,e,work,a,info)
      return
      end
