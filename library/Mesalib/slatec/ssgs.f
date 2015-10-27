*deck ssgs
      subroutine ssgs (n, b, x, nelt, ia, ja, a, isym, itol, tol, itmax,
     +   iter, err, ierr, iunit, rwork, lenw, iwork, leniw)
c***begin prologue  ssgs
c***purpose  gauss-seidel method iterative sparse ax = b solver.
c            routine to solve a general linear system  ax = b  using
c            gauss-seidel iteration.
c***library   slatec (slap)
c***category  d2a4, d2b4
c***type      single precision (ssgs-s, dsgs-d)
c***keywords  iterative precondition, linear system, slap, sparse
c***author  greenbaum, anne, (courant institute)
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c *usage:
c     integer n, nelt, ia(nelt), ja(nelt), isym, itol, itmax
c     integer iter, ierr, iunit, lenw, iwork(nl+2*n+1), leniw
c     real b(n), x(n), a(nelt), tol, err, rwork(nl+3*n)
c
c     call ssgs(n, b, x, nelt, ia, ja, a, isym, itol, tol,
c    $     itmax, iter, err, ierr, iunit, rwork, lenw, iwork, leniw )
c
c *arguments:
c n      :in       integer.
c         order of the matrix.
c b      :in       real b(n).
c         right-hand side vector.
c x      :inout    real x(n).
c         on input x is your initial guess for solution vector.
c         on output x is the final approximate solution.
c nelt   :in       integer.
c         number of non-zeros stored in a.
c ia     :inout    integer ia(nelt).
c ja     :inout    integer ja(nelt).
c a      :inout    real a(nelt).
c         these arrays should hold the matrix a in either the slap
c         triad format or the slap column format.  see "description",
c         below.  if the slap triad format is chosen it is changed
c         internally to the slap column format.
c isym   :in       integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all non-zero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the lower
c         lower triangle of the matrix is stored.
c itol   :in       integer.
c         flag to indicate type of convergence criterion.
c         if itol=1, iteration stops when the 2-norm of the residual
c         divided by the 2-norm of the right-hand side is less than tol.
c         if itol=2, iteration stops when the 2-norm of m-inv times the
c         residual divided by the 2-norm of m-inv times the right hand
c         side is less than tol, where m-inv is the inverse of the
c         diagonal of a.
c         itol=11 is often useful for checking and comparing different
c         routines.  for this case, the user must supply the "exact"
c         solution or a very accurate approximation (one with an error
c         much less than tol) through a common block,
c             common /sslblk/ soln( )
c         if itol=11, iteration stops when the 2-norm of the difference
c         between the iterative approximation and the user-supplied
c         solution divided by the 2-norm of the user-supplied solution
c         is less than tol.  note that this requires the user to set up
c         the "common /sslblk/ soln(length)" in the calling routine.
c         the routine with this declaration should be loaded before the
c         stop test so that the correct length is used by the loader.
c         this procedure is not standard fortran and may not work
c         correctly on your system (although it has worked on every
c         system the authors have tried).  if itol is not 11 then this
c         common block is indeed standard fortran.
c tol    :inout    real.
c         convergence criterion, as described above.  (reset if ierr=4.)
c itmax  :in       integer.
c         maximum number of iterations.
c iter   :out      integer.
c         number of iterations required to reach convergence, or
c         itmax+1 if convergence criterion could not be achieved in
c         itmax iterations.
c err    :out      real.
c         error estimate of error in final approximate solution, as
c         defined by itol.
c ierr   :out      integer.
c         return error flag.
c           ierr = 0 => all went well.
c           ierr = 1 => insufficient space allocated for work or iwork.
c           ierr = 2 => method failed to converge in itmax steps.
c           ierr = 3 => error in user input.
c                       check input values of n, itol.
c           ierr = 4 => user error tolerance set too tight.
c                       reset to 500*r1mach(3).  iteration proceeded.
c           ierr = 5 => preconditioning matrix, m, is not positive
c                       definite.  (r,z) < 0.
c           ierr = 6 => matrix a is not positive definite.  (p,ap) < 0.
c iunit  :in       integer.
c         unit number on which to write the error at each iteration,
c         if this is desired for monitoring convergence.  if unit
c         number is 0, no writing will occur.
c rwork  :work     real rwork(lenw).
c         real array used for workspace.
c lenw   :in       integer.
c         length of the real workspace, rwork.  lenw >= nl+3*n.
c         nl is the number of non-zeros in the lower triangle of the
c         matrix (including the diagonal).
c iwork  :work     integer iwork(leniw).
c         integer array used for workspace.
c         upon return the following locations of iwork hold information
c         which may be of use to the user:
c         iwork(9)  amount of integer workspace actually used.
c         iwork(10) amount of real workspace actually used.
c leniw  :in       integer.
c         length of the integer workspace, iwork.  leniw >= nl+n+11.
c         nl is the number of non-zeros in the lower triangle of the
c         matrix (including the diagonal).
c
c *description
c       the sparse linear algebra package (slap) utilizes two matrix
c       data structures: 1) the  slap triad  format or  2)  the slap
c       column format.  the user can hand this routine either of the
c       of these data structures and slap  will figure out  which on
c       is being used and act accordingly.
c
c       =================== s l a p triad format ===================
c
c       this routine requires that the  matrix a be   stored in  the
c       slap  triad format.  in  this format only the non-zeros  are
c       stored.  they may appear in  *any* order.  the user supplies
c       three arrays of  length nelt, where  nelt is  the number  of
c       non-zeros in the matrix: (ia(nelt), ja(nelt), a(nelt)).  for
c       each non-zero the user puts the row and column index of that
c       matrix element  in the ia and  ja arrays.  the  value of the
c       non-zero   matrix  element is  placed  in  the corresponding
c       location of the a array.   this is  an  extremely  easy data
c       structure to generate.  on  the  other hand it   is  not too
c       efficient on vector computers for  the iterative solution of
c       linear systems.  hence,   slap changes   this  input    data
c       structure to the slap column format  for  the iteration (but
c       does not change it back).
c
c       here is an example of the  slap triad   storage format for a
c       5x5 matrix.  recall that the entries may appear in any order.
c
c           5x5 matrix      slap triad format for 5x5 matrix on left.
c                              1  2  3  4  5  6  7  8  9 10 11
c       |11 12  0  0 15|   a: 51 12 11 33 15 53 55 22 35 44 21
c       |21 22  0  0  0|  ia:  5  1  1  3  1  5  5  2  3  4  2
c       | 0  0 33  0 35|  ja:  1  2  1  3  5  3  5  2  5  4  1
c       | 0  0  0 44  0|
c       |51  0 53  0 55|
c
c       =================== s l a p column format ==================
c
c       this routine  requires that  the matrix a  be stored in  the
c       slap column format.  in this format the non-zeros are stored
c       counting down columns (except for  the diagonal entry, which
c       must appear first in each  "column")  and are stored  in the
c       real array a.  in other words, for each column in the matrix
c       put the diagonal entry in a.  then put in the other non-zero
c       elements going down   the  column (except  the diagonal)  in
c       order.  the ia array holds the row  index for each non-zero.
c       the ja array holds the offsets into the ia, a arrays for the
c       beginning of   each    column.    that  is,    ia(ja(icol)),
c       a(ja(icol)) points to the beginning of the icol-th column in
c       ia and  a.  ia(ja(icol+1)-1),  a(ja(icol+1)-1) points to the
c       end  of   the icol-th  column.  note   that  we  always have
c       ja(n+1) = nelt+1, where  n  is the number of columns in  the
c       matrix and  nelt   is the number of non-zeros in the matrix.
c
c       here is an example of the  slap column  storage format for a
c       5x5 matrix (in the a and ia arrays '|'  denotes the end of a
c       column):
c
c           5x5 matrix      slap column format for 5x5 matrix on left.
c                              1  2  3    4  5    6  7    8    9 10 11
c       |11 12  0  0 15|   a: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
c       |21 22  0  0  0|  ia:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
c       | 0  0 33  0 35|  ja:  1  4  6    8  9   12
c       | 0  0  0 44  0|
c       |51  0 53  0 55|
c
c *side effects:
c       the slap triad format (ia, ja, a) is modified internally to be
c       the slap column format.  see above.
c
c *cautions:
c     this routine will attempt to write to the fortran logical output
c     unit iunit, if iunit .ne. 0.  thus, the user must make sure that
c     this logical unit is attached to a file or terminal before calling
c     this routine with a non-zero value for iunit.  this routine does
c     not check for the validity of a non-zero iunit unit number.
c
c***see also  ssjac, sir
c***references  (none)
c***routines called  schkw, sir, ss2lt, ss2y, ssli, ssmv
c***revision history  (yymmdd)
c   871119  date written
c   881213  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890921  removed tex from comments.  (fnf)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   910411  prologue converted to version 4.0 format.  (bab)
c   920407  common block renamed sslblk.  (wrb)
c   920511  added complete declaration section.  (wrb)
c   921019  corrected nel to nl.  (fnf)
c***end prologue  ssgs
c     .. parameters ..
      integer locrb, locib
      parameter (locrb=1, locib=11)
c     .. scalar arguments ..
      real err, tol
      integer ierr, isym, iter, itmax, itol, iunit, leniw, lenw, n, nelt
c     .. array arguments ..
      real a(n), b(n), rwork(*), x(n)
      integer ia(nelt), iwork(10), ja(nelt)
c     .. local scalars ..
      integer icol, j, jbgn, jend, locdz, locel, lociel, lociw, locjel,
     +        locr, locw, locz, nl
c     .. external subroutines ..
      external schkw, sir, ss2lt, ss2y, ssli, ssmv
c***first executable statement  ssgs
c
      if( n.lt.1 .or. nelt.lt.1 ) then
         ierr = 3
         return
      endif
c
c         modify the slap matrix data structure to ysmp-column.
      call ss2y( n, nelt, ia, ja, a, isym )
c
c         count number of elements in lower triangle of the matrix.
      if( isym.eq.0 ) then
         nl = 0
         do 20 icol = 1, n
            jbgn = ja(icol)
            jend = ja(icol+1)-1
            do 10 j = jbgn, jend
               if( ia(j).ge.icol ) nl = nl + 1
 10         continue
 20      continue
      else
         nl = ja(n+1)-1
      endif
c
c         set up the work arrays.  then store the lower triangle of
c         the matrix.
c
      locjel = locib
      lociel = locjel + n+1
      lociw = lociel + nl
c
      locel = locrb
      locr = locel + nl
      locz = locr + n
      locdz = locz + n
      locw = locdz + n
c
c         check the workspace allocations.
      call schkw( 'ssgs', lociw, leniw, locw, lenw, ierr, iter, err )
      if( ierr.ne.0 ) return
c
      iwork(1) = nl
      iwork(2) = lociel
      iwork(3) = locjel
      iwork(4) = locel
      iwork(9) = lociw
      iwork(10) = locw
c
      call ss2lt( n, nelt, ia, ja, a, isym, nl, iwork(lociel),
     $     iwork(locjel), rwork(locel) )
c
c         call iterative refinement routine.
      call sir(n, b, x, nelt, ia, ja, a, isym, ssmv, ssli,
     $     itol, tol, itmax, iter, err, ierr, iunit, rwork(locr),
     $     rwork(locz), rwork(locdz), rwork, iwork )
c
c         set the amount of integer and real workspace used.
      iwork(9) = lociw+n+nelt
      iwork(10) = locw+nelt
      return
c------------- last line of ssgs follows ------------------------------
      end
