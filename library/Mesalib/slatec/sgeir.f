*deck sgeir
      subroutine sgeir (a, lda, n, v, itask, ind, work, iwork)
c***begin prologue  sgeir
c***purpose  solve a general system of linear equations.  iterative
c            refinement is used to obtain an error estimate.
c***library   slatec
c***category  d2a1
c***type      single precision (sgeir-s, cgeir-c)
c***keywords  complex linear equations, general matrix,
c             general system of linear equations
c***author  voorhees, e. a., (lanl)
c***description
c
c    subroutine sgeir solves a general nxn system of single
c    precision linear equations using linpack subroutines sgefa and
c    sgesl.  one pass of iterative refinement is used only to obtain
c    an estimate of the accuracy.  that is, if a is an nxn real
c    matrix and if x and b are real n-vectors, then sgeir solves
c    the equation
c
c                          a*x=b.
c
c    the matrix a is first factored into upper and lower tri-
c    angular matrices u and l using partial pivoting.  these
c    factors and the pivoting information are used to calculate
c    the solution, x.  then the residual vector is found and
c    used to calculate an estimate of the relative error, ind.
c    ind estimates the accuracy of the solution only when the
c    input matrix and the right hand side are represented
c    exactly in the computer and does not take into account
c    any errors in the input data.
c
c    if the equation a*x=b is to be solved for more than one vector
c    b, the factoring of a does not need to be performed again and
c    the option to solve only (itask .gt. 1) will be faster for
c    the succeeding solutions.  in this case, the contents of a,
c    lda, n, work, and iwork must not have been altered by the
c    user following factorization (itask=1).  ind will not be
c    changed by sgeir in this case.
c
c  argument description ***
c
c    a      real(lda,n)
c             the doubly subscripted array with dimension (lda,n)
c             which contains the coefficient matrix.  a is not
c             altered by the routine.
c    lda    integer
c             the leading dimension of the array a.  lda must be great-
c             er than or equal to n.  (terminal error message ind=-1)
c    n      integer
c             the order of the matrix a.  the first n elements of
c             the array a are the elements of the first column of
c             matrix a.  n must be greater than or equal to 1.
c             (terminal error message ind=-2)
c    v      real(n)
c             on entry, the singly subscripted array(vector) of di-
c               mension n which contains the right hand side b of a
c               system of simultaneous linear equations a*x=b.
c             on return, v contains the solution vector, x .
c    itask  integer
c             if itask=1, the matrix a is factored and then the
c               linear equation is solved.
c             if itask .gt. 1, the equation is solved using the existing
c               factored matrix a (stored in work).
c             if itask .lt. 1, then terminal error message ind=-3 is
c               printed.
c    ind    integer
c             gt. 0  ind is a rough estimate of the number of digits
c                     of accuracy in the solution, x.  ind=75 means
c                     that the solution vector x is zero.
c             lt. 0  see error message corresponding to ind below.
c    work   real(n*(n+1))
c             a singly subscripted array of dimension at least n*(n+1).
c    iwork  integer(n)
c             a singly subscripted array of dimension at least n.
c
c  error messages printed ***
c
c    ind=-1  terminal   n is greater than lda.
c    ind=-2  terminal   n is less than one.
c    ind=-3  terminal   itask is less than one.
c    ind=-4  terminal   the matrix a is computationally singular.
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
c***routines called  r1mach, sasum, scopy, sdsdot, sgefa, sgesl, xermsg
c***revision history  (yymmdd)
c   800430  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  sgeir
c
      integer lda,n,itask,ind,iwork(*),info,j
      real a(lda,*),v(*),work(n,*),xnorm,dnorm,sdsdot,sasum,r1mach
      character*8 xern1, xern2
c***first executable statement  sgeir
      if (lda.lt.n) then
         ind = -1
         write (xern1, '(i8)') lda
         write (xern2, '(i8)') n
         call xermsg ('slatec', 'sgeir', 'lda = ' // xern1 //
     *      ' is less than n = ' // xern2, -1, 1)
         return
      endif
c
      if (n.le.0) then
         ind = -2
         write (xern1, '(i8)') n
         call xermsg ('slatec', 'sgeir', 'n = ' // xern1 //
     *      ' is less than 1', -2, 1)
         return
      endif
c
      if (itask.lt.1) then
         ind = -3
         write (xern1, '(i8)') itask
         call xermsg ('slatec', 'sgeir', 'itask = ' // xern1 //
     *      ' is less than 1', -3, 1)
         return
      endif
c
      if (itask.eq.1) then
c
c        move matrix a to work
c
         do 10 j=1,n
            call scopy(n,a(1,j),1,work(1,j),1)
   10    continue
c
c        factor matrix a into lu
c
         call sgefa(work,n,n,iwork,info)
c
c        check for computationally singular matrix
c
         if (info.ne.0) then
            ind = -4
            call xermsg ('slatec', 'sgeir',
     *         'singular matrix a - no solution', -4, 1)
            return
         endif
      endif
c
c     solve when factoring complete
c     move vector b to work
c
      call scopy(n,v(1),1,work(1,n+1),1)
      call sgesl(work,n,n,iwork,v,0)
c
c     form norm of x0
c
      xnorm=sasum(n,v(1),1)
      if (xnorm.ne.0.0) then
         ind = 75
         return
      endif
c
c     compute  residual
c
      do 40 j=1,n
         work(j,n+1) = sdsdot(n,-work(j,n+1),a(j,1),lda,v,1)
   40 continue
c
c     solve a*delta=r
c
      call sgesl(work,n,n,iwork,work(1,n+1),0)
c
c     form norm of delta
c
      dnorm = sasum(n,work(1,n+1),1)
c
c     compute ind (estimate of no. of significant digits)
c     and check for ind greater than zero
c
      ind = -log10(max(r1mach(4),dnorm/xnorm))
      if (ind.le.0) then
         ind = -10
         call xermsg ('slatec', 'sgeir',
     *      'solution may have no significance', -10, 0)
      endif
      return
      end
