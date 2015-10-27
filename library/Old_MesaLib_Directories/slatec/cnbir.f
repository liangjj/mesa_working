*deck cnbir
      subroutine cnbir (abe, lda, n, ml, mu, v, itask, ind, work, iwork)
c***begin prologue  cnbir
c***purpose  solve a general nonsymmetric banded system of linear
c            equations.  iterative refinement is used to obtain an error
c            estimate.
c***library   slatec
c***category  d2c2
c***type      complex (snbir-s, cnbir-c)
c***keywords  banded, linear equations, nonsymmetric
c***author  voorhees, e. a., (lanl)
c***description
c
c    subroutine cnbir solves a general nonsymmetric banded nxn
c    system of single precision complex linear equations using
c    slatec subroutines cnbfa and cnbsl.  these are adaptations
c    of the linpack subroutines cgbfa and cgbsl which require
c    a different format for storing the matrix elements.
c    one pass of iterative refinement is used only to obtain an
c    estimate of the accuracy.  if  a  is an nxn complex banded
c    matrix and if  x  and  b  are complex n-vectors, then cnbir
c    solves the equation
c
c                          a*x=b.
c
c    a band matrix is a matrix whose nonzero elements are all
c    fairly near the main diagonal, specifically  a(i,j) = 0
c    if  i-j is greater than  ml  or  j-i  is greater than
c    mu .  the integers ml and mu are called the lower and upper
c    band widths and  m = ml+mu+1  is the total band width.
c    cnbir uses less time and storage than the corresponding
c    program for general matrices (cgeir) if 2*ml+mu .lt. n .
c
c    the matrix a is first factored into upper and lower tri-
c    angular matrices u and l using partial pivoting.  these
c    factors and the pivoting information are used to find the
c    solution vector x .  then the residual vector is found and used
c    to calculate an estimate of the relative error, ind .  ind esti-
c    mates the accuracy of the solution only when the input matrix
c    and the right hand side are represented exactly in the computer
c    and does not take into account any errors in the input data.
c
c    if the equation a*x=b is to be solved for more than one vector
c    b, the factoring of a does not need to be performed again and
c    the option to only solve (itask .gt. 1) will be faster for
c    the succeeding solutions.  in this case, the contents of a, lda,
c    n, work and iwork must not have been altered by the user follow-
c    ing factorization (itask=1).  ind will not be changed by cnbir
c    in this case.
c
c
c    band storage
c
c          if  a  is a band matrix, the following program segment
c          will set up the input.
c
c                  ml = (band width below the diagonal)
c                  mu = (band width above the diagonal)
c                  do 20 i = 1, n
c                     j1 = max(1, i-ml)
c                     j2 = min(n, i+mu)
c                     do 10 j = j1, j2
c                        k = j - i + ml + 1
c                        abe(i,k) = a(i,j)
c               10    continue
c               20 continue
c
c          this uses columns  1  through  ml+mu+1  of abe .
c
c    example:  if the original matrix is
c
c          11 12 13  0  0  0
c          21 22 23 24  0  0
c           0 32 33 34 35  0
c           0  0 43 44 45 46
c           0  0  0 54 55 56
c           0  0  0  0 65 66
c
c     then  n = 6, ml = 1, mu = 2, lda .ge. 5  and abe should contain
c
c           * 11 12 13        , * = not used
c          21 22 23 24
c          32 33 34 35
c          43 44 45 46
c          54 55 56  *
c          65 66  *  *
c
c
c  argument description ***
c
c    abe    complex(lda,mm)
c             on entry, contains the matrix in band storage as
c               described above.  mm  must not be less than  m =
c               ml+mu+1 .  the user is cautioned to dimension  abe
c               with care since mm is not an argument and cannot
c               be checked by cnbir.  the rows of the original
c               matrix are stored in the rows of  abe  and the
c               diagonals of the original matrix are stored in
c               columns  1  through  ml+mu+1  of  abe .  abe  is
c               not altered by the program.
c    lda    integer
c             the leading dimension of array abe.  lda must be great-
c             er than or equal to n.  (terminal error message ind=-1)
c    n      integer
c             the order of the matrix a.  n must be greater
c             than or equal to 1 .  (terminal error message ind=-2)
c    ml     integer
c             the number of diagonals below the main diagonal.
c             ml  must not be less than zero nor greater than or
c             equal to  n .  (terminal error message ind=-5)
c    mu     integer
c             the number of diagonals above the main diagonal.
c             mu  must not be less than zero nor greater than or
c             equal to  n .  (terminal error message ind=-6)
c    v      complex(n)
c             on entry, the singly subscripted array(vector) of di-
c               mension n which contains the right hand side b of a
c               system of simultaneous linear equations a*x=b.
c             on return, v contains the solution vector, x .
c    itask  integer
c             if itask=1, the matrix a is factored and then the
c               linear equation is solved.
c             if itask .gt. 1, the equation is solved using the existing
c               factored matrix a and iwork.
c             if itask .lt. 1, then terminal error message ind=-3 is
c               printed.
c    ind    integer
c             gt. 0  ind is a rough estimate of the number of digits
c                     of accuracy in the solution, x .  ind=75 means
c                     that the solution vector  x  is zero.
c             lt. 0  see error message corresponding to ind below.
c    work   complex(n*(nc+1))
c             a singly subscripted array of dimension at least
c             n*(nc+1)  where  nc = 2*ml+mu+1 .
c    iwork  integer(n)
c             a singly subscripted array of dimension at least n.
c
c  error messages printed ***
c
c    ind=-1  terminal   n is greater than lda.
c    ind=-2  terminal   n is less than 1.
c    ind=-3  terminal   itask is less than 1.
c    ind=-4  terminal   the matrix a is computationally singular.
c                         a solution has not been computed.
c    ind=-5  terminal   ml is less than zero or is greater than
c                         or equal to n .
c    ind=-6  terminal   mu is less than zero or is greater than
c                         or equal to n .
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
c***routines called  ccopy, cdcdot, cnbfa, cnbsl, r1mach, scasum, xermsg
c***revision history  (yymmdd)
c   800819  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900510  convert xerrwv calls to xermsg calls, cvt goto's to
c           if-then-else.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cnbir
c
      integer lda,n,itask,ind,iwork(*),info,j,k,kk,l,m,ml,mu,nc
      complex abe(lda,*),v(*),work(n,*),cdcdot
      real xnorm,dnorm,scasum,r1mach
      character*8 xern1, xern2
c***first executable statement  cnbir
      if (lda.lt.n) then
         ind = -1
         write (xern1, '(i8)') lda
         write (xern2, '(i8)') n
         call xermsg ('slatec', 'cnbir', 'lda = ' // xern1 //
     *      ' is less than n = ' // xern2, -1, 1)
         return
      endif
c
      if (n.le.0) then
         ind = -2
         write (xern1, '(i8)') n
         call xermsg ('slatec', 'cnbir', 'n = ' // xern1 //
     *      ' is less than 1', -2, 1)
         return
      endif
c
      if (itask.lt.1) then
         ind = -3
         write (xern1, '(i8)') itask
         call xermsg ('slatec', 'cnbir', 'itask = ' // xern1 //
     *      ' is less than 1', -3, 1)
         return
      endif
c
      if (ml.lt.0 .or. ml.ge.n) then
         ind = -5
         write (xern1, '(i8)') ml
         call xermsg ('slatec', 'cnbir',
     *      'ml = ' // xern1 // ' is out of range', -5, 1)
         return
      endif
c
      if (mu.lt.0 .or. mu.ge.n) then
         ind = -6
         write (xern1, '(i8)') mu
         call xermsg ('slatec', 'cnbir',
     *      'mu = ' // xern1 // ' is out of range', -6, 1)
         return
      endif
c
      nc = 2*ml+mu+1
      if (itask.eq.1) then
c
c        move matrix abe to work
c
         m=ml+mu+1
         do 10 j=1,m
            call ccopy(n,abe(1,j),1,work(1,j),1)
   10    continue
c
c        factor matrix a into lu
         call cnbfa(work,n,n,ml,mu,iwork,info)
c
c        check for computationally singular matrix
         if (info.ne.0) then
            ind=-4
            call xermsg ('slatec', 'cnbir',
     *         'singular matrix a - no solution', -4, 1)
            return
         endif
      endif
c
c     solve when factoring complete
c     move vector b to work
c
      call ccopy(n,v(1),1,work(1,nc+1),1)
      call cnbsl(work,n,n,ml,mu,iwork,v,0)
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
         k  = max(1,ml+2-j)
         kk = max(1,j-ml)
         l  = min(j-1,ml)+min(n-j,mu)+1
         work(j,nc+1) = cdcdot(l,-work(j,nc+1),abe(j,k),lda,v(kk),1)
   40 continue
c
c     solve a*delta=r
c
      call cnbsl(work,n,n,ml,mu,iwork,work(1,nc+1),0)
c
c     form norm of delta
c
      dnorm = scasum(n,work(1,nc+1),1)
c
c     compute ind (estimate of no. of significant digits)
c     and check for ind greater than zero
c
      ind = -log10(max(r1mach(4),dnorm/xnorm))
      if (ind.le.0) then
         ind = -10
         call xermsg ('slatec', 'cnbir',
     *      'solution may have no significance', -10, 0)
      endif
      return
      end
