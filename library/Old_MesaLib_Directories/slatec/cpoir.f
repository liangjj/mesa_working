*deck cpoir
      subroutine cpoir (a, lda, n, v, itask, ind, work)
c***begin prologue  cpoir
c***purpose  solve a positive definite hermitian system of linear
c            equations.  iterative refinement is used to obtain an
c            error estimate.
c***library   slatec
c***category  d2d1b
c***type      complex (spoir-s, cpoir-c)
c***keywords  hermitian, linear equations, positive definite, symmetric
c***author  voorhees, e. a., (lanl)
c***description
c
c    subroutine cpoir solves a complex positive definite hermitian
c    nxn system of single precision linear equations using linpack
c    subroutines cpofa and cposl.  one pass of iterative refine-
c    ment is used only to obtain an estimate of the accuracy.  that
c    is, if a is an nxn complex positive definite hermitian matrix
c    and if x and b are complex n-vectors, then cpoir solves the
c    equation
c
c                          a*x=b.
c
c    care should be taken not to use cpoir with a non-hermitian
c    matrix.
c
c    the matrix a is first factored into upper and lower
c    triangular matrices r and r-transpose.  these
c    factors are used to calculate the solution, x.
c    then the residual vector is found and used
c    to calculate an estimate of the relative error, ind.
c    ind estimates the accuracy of the solution only when the
c    input matrix and the right hand side are represented
c    exactly in the computer and does not take into account
c    any errors in the input data.
c
c    if the equation a*x=b is to be solved for more than one vector
c    b, the factoring of a does not need to be performed again and
c    the option to only solve (itask .gt. 1) will be faster for
c    the succeeding solutions.  in this case, the contents of a,
c    lda, n, and work must not have been altered by the user
c    following factorization (itask=1).  ind will not be changed
c    by cpoir in this case.
c
c  argument description ***
c    a       complex(lda,n)
c             the doubly subscripted array with dimension (lda,n)
c             which contains the coefficient matrix.  only the
c             upper triangle, including the diagonal, of the
c             coefficient matrix need be entered.  a is not
c             altered by the routine.
c    lda    integer
c             the leading dimension of the array a.  lda must be great-
c             er than or equal to n.  (terminal error message ind=-1)
c    n      integer
c             the order of the matrix a.  n must be greater than
c             or equal to one.  (terminal error message ind=-2)
c    v      complex(n)
c             on entry, the singly subscripted array(vector) of di-
c               mension n which contains the right hand side b of a
c               system of simultaneous linear equations a*x=b.
c             on return, v contains the solution vector, x .
c    itask  integer
c             if itask = 1, the matrix a is factored and then the
c               linear equation is solved.
c             if itask .gt. 1, the equation is solved using the existing
c               factored matrix a (stored in work).
c             if itask .lt. 1, then terminal terminal error ind=-3 is
c               printed.
c    ind    integer
c             gt. 0  ind is a rough estimate of the number of digits
c                     of accuracy in the solution, x.  ind=75 means
c                     that the solution vector x is zero.
c             lt. 0  see error message corresponding to ind below.
c    work   complex(n*(n+1))
c             a singly subscripted array of dimension at least n*(n+1).
c
c  error messages printed ***
c
c    ind=-1  terminal   n is greater than lda.
c    ind=-2  terminal   n is less than one.
c    ind=-3  terminal   itask is less than one.
c    ind=-4  terminal   the matrix a is computationally singular
c                         or is not positive definite.
c                         a solution has not been computed.
c    ind=-10 warning    the solution has no apparent significance.
c                         the solution may be inaccurate or the matrix
c                         a may be poorly scaled.
c
c               note-  the above terminal(*fatal*) error messages are
c                      designed to be handled by xermsg in which
c                      level=1 (recoverable) and iflag=2 .  level=0
c                      for warning error messages from xermsg.  unless
c                      the user provides otherwise, an error message
c                      will be printed followed by an abort.
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  ccopy, cpofa, cposl, dcdot, r1mach, scasum, xermsg
c***revision history  (yymmdd)
c   800530  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900510  convert xerrwv calls to xermsg calls, cvt goto's to
c           if-then-else.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cpoir
c
      integer lda,n,itask,ind,info,j
      complex a(lda,*),v(*),work(n,*)
      real scasum,xnorm,dnorm,r1mach
      double precision dr1,di1,dr2,di2
      character*8 xern1, xern2
c***first executable statement  cpoir
      if (lda.lt.n) then
         ind = -1
         write (xern1, '(i8)') lda
         write (xern2, '(i8)') n
         call xermsg ('slatec', 'cpoir', 'lda = ' // xern1 //
     *      ' is less than n = ' // xern2, -1, 1)
         return
      endif
c
      if (n.le.0) then
         ind = -2
         write (xern1, '(i8)') n
         call xermsg ('slatec', 'cpoir', 'n = ' // xern1 //
     *      ' is less than 1', -2, 1)
         return
      endif
c
      if (itask.lt.1) then
         ind = -3
         write (xern1, '(i8)') itask
         call xermsg ('slatec', 'cpoir', 'itask = ' // xern1 //
     *      ' is less than 1', -3, 1)
         return
      endif
c
      if (itask.eq.1) then
c
c        move matrix a to work
c
         do 10 j=1,n
            call ccopy(n,a(1,j),1,work(1,j),1)
   10    continue
c
c        factor matrix a into r
c
         call cpofa(work,n,n,info)
c
c        check for  singular or not pos.def. matrix
c
         if (info.ne.0) then
            ind = -4
            call xermsg ('slatec', 'cpoir',
     *         'singular or not positive definite - no solution', -4, 1)
            return
         endif
      endif
c
c     solve after factoring
c     move vector b to work
c
      call ccopy(n,v(1),1,work(1,n+1),1)
      call cposl(work,n,n,v)
c
c     form norm of x0
c
      xnorm = scasum(n,v(1),1)
      if (xnorm.eq.0.0) then
         ind = 75
         return
      endif
c
c     compute  residual
c
      do 40 j=1,n
         call dcdot(j-1,-1.d0,a(1,j),1,v(1),1,dr1,di1)
         call dcdot(n-j+1,1.d0,a(j,j),lda,v(j),1,dr2,di2)
         dr1 = dr1+dr2-dble(real(work(j,n+1)))
         di1 = di1+di2-dble(aimag(work(j,n+1)))
         work(j,n+1) = cmplx(real(dr1),real(di1))
   40 continue
c
c     solve a*delta=r
c
      call cposl(work,n,n,work(1,n+1))
c
c     form norm of delta
c
      dnorm = scasum(n,work(1,n+1),1)
c
c     compute ind (estimate of no. of significant digits)
c     and check for ind greater than zero
c
      ind = -log10(max(r1mach(4),dnorm/xnorm))
      if (ind.le.0) then
         ind = -10
         call xermsg ('slatec', 'cpoir',
     *      'solution may have no significance', -10, 0)
      endif
      return
      end
