*deck dsdcg
      subroutine dsdcg (n, b, x, nelt, ia, ja, a, isym, itol, tol,
     +   itmax, iter, err, ierr, iunit, rwork, lenw, iwork, leniw)
c***begin prologue  dsdcg
c***purpose  diagonally scaled conjugate gradient sparse ax=b solver.
c            routine to solve a symmetric positive definite linear
c            system  ax = b  using the preconditioned conjugate
c            gradient method.  the preconditioner is diagonal scaling.
c***library   slatec (slap)
c***category  d2b4
c***type      double precision (ssdcg-s, dsdcg-d)
c***keywords  iterative precondition, slap, sparse,
c             symmetric linear system
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
c     integer iter, ierr, iunit, lenw, iwork(10), leniw
c     double precision b(n), x(n), a(nelt), tol, err, rwork(5*n)
c
c     call dsdcg(n, b, x, nelt, ia, ja, a, isym, itol, tol,
c    $     itmax, iter, err, ierr, iunit, rwork, lenw, iwork, leniw )
c
c *arguments:
c n      :in       integer.
c         order of the matrix.
c b      :in       double precision b(n).
c         right-hand side vector.
c x      :inout    double precision x(n).
c         on input x is your initial guess for solution vector.
c         on output x is the final approximate solution.
c nelt   :in       integer.
c         number of non-zeros stored in a.
c ia     :inout    integer ia(nelt).
c ja     :inout    integer ja(nelt).
c a      :inout    double precision a(nelt).
c         these arrays should hold the matrix a in either the slap
c         triad format or the slap column format.  see "description",
c         below.  if the slap triad format is chosen it is changed
c         internally to the slap column format.
c isym   :in       integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all non-zero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the upper
c         or lower triangle of the matrix is stored.
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
c             common /dslblk/ soln( )
c         if itol=11, iteration stops when the 2-norm of the difference
c         between the iterative approximation and the user-supplied
c         solution divided by the 2-norm of the user-supplied solution
c         is less than tol.  note that this requires the user to set up
c         the "common /dslblk/ soln(length)" in the calling routine.
c         the routine with this declaration should be loaded before the
c         stop test so that the correct length is used by the loader.
c         this procedure is not standard fortran and may not work
c         correctly on your system (although it has worked on every
c         system the authors have tried).  if itol is not 11 then this
c         common block is indeed standard fortran.
c tol    :inout    double precision.
c         convergence criterion, as described above.  (reset if ierr=4.)
c itmax  :in       integer.
c         maximum number of iterations.
c iter   :out      integer.
c         number of iterations required to reach convergence, or
c         itmax+1 if convergence criterion could not be achieved in
c         itmax iterations.
c err    :out      double precision.
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
c                       reset to 500*d1mach(3).  iteration proceeded.
c           ierr = 5 => preconditioning matrix, m, is not positive
c                       definite.  (r,z) < 0.
c           ierr = 6 => matrix a is not positive definite.  (p,ap) < 0.
c iunit  :in       integer.
c         unit number on which to write the error at each iteration,
c         if this is desired for monitoring convergence.  if unit
c         number is 0, no writing will occur.
c rwork  :work     double precision rwork(lenw).
c         double precision array used for workspace.
c lenw   :in       integer.
c         length of the double precision workspace, rwork.  lenw >= 5*n.
c iwork  :work     integer iwork(leniw).
c         used to hold pointers into the double precision workspace,
c         rwork.  upon return the following locations of iwork hold
c         information which may be of use to the user:
c         iwork(9)  amount of integer workspace actually used.
c         iwork(10) amount of double precision workspace actually used.
c leniw  :in       integer.
c         length of the integer workspace, iwork.  leniw >= 10.
c
c *description:
c       this  routine   performs preconditioned conjugate   gradient
c       method on  the  symmetric positive definite   linear  system
c       ax=b.   the preconditioner is  m = diag(a), the  diagonal of
c       the matrix a.  this is the  simplest of preconditioners  and
c       vectorizes very well.   this routine is  simply a driver for
c       the dcg routine.  it  calls the dsds  routine to  set up the
c       preconditioning  and  then  calls  dcg  with the appropriate
c       matvec and msolve routines.
c
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
c       double precision array a.   in other words,  for each column
c       in the matrix put the diagonal entry in  a.  then put in the
c       other non-zero  elements going down  the column (except  the
c       diagonal) in order.   the  ia array holds the  row index for
c       each non-zero.  the ja array holds the offsets  into the ia,
c       a arrays  for  the  beginning  of each   column.   that  is,
c       ia(ja(icol)),  a(ja(icol)) points   to the beginning  of the
c       icol-th   column    in    ia and   a.      ia(ja(icol+1)-1),
c       a(ja(icol+1)-1) points to  the  end of the   icol-th column.
c       note that we always have  ja(n+1) = nelt+1,  where n is  the
c       number of columns in  the matrix and nelt  is the number  of
c       non-zeros in the matrix.
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
c       the slap triad format (ia, ja, a) is modified internally  to
c       be the slap column format.  see above.
c
c *cautions:
c     this routine will attempt to write to the fortran logical output
c     unit iunit, if iunit .ne. 0.  thus, the user must make sure that
c     this logical unit is attached to a file or terminal before calling
c     this routine with a non-zero value for iunit.  this routine does
c     not check for the validity of a non-zero iunit unit number.
c
c***see also  dcg, dsiccg
c***references  1. louis hageman and david young, applied iterative
c                  methods, academic press, new york, 1981.
c               2. concus, golub and o'leary, a generalized conjugate
c                  gradient method for the numerical solution of
c                  elliptic partial differential equations, in sparse
c                  matrix computations, bunch and rose, eds., academic
c                  press, new york, 1979.
c***routines called  dcg, dchkw, ds2y, dsdi, dsds, dsmv
c***revision history  (yymmdd)
c   890404  date written
c   890404  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890921  removed tex from comments.  (fnf)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   910411  prologue converted to version 4.0 format.  (bab)
c   920407  common block renamed dslblk.  (wrb)
c   920511  added complete declaration section.  (wrb)
c   920929  corrected format of references.  (fnf)
c***end prologue  dsdcg
c     .. parameters ..
      integer locrb, locib
      parameter (locrb=1, locib=11)
c     .. scalar arguments ..
      double precision err, tol
      integer ierr, isym, iter, itmax, itol, iunit, leniw, lenw, n, nelt
c     .. array arguments ..
      double precision a(nelt), b(n), rwork(lenw), x(n)
      integer ia(nelt), iwork(leniw), ja(nelt)
c     .. local scalars ..
      integer locd, locdz, lociw, locp, locr, locw, locz
c     .. external subroutines ..
      external dcg, dchkw, ds2y, dsdi, dsds, dsmv
c***first executable statement  dsdcg
c
      ierr = 0
      if( n.lt.1 .or. nelt.lt.1 ) then
         ierr = 3
         return
      endif
c
c         modify the slap matrix data structure to ysmp-column.
      call ds2y( n, nelt, ia, ja, a, isym )
c
c         set up the work arrays.
      lociw = locib
c
      locd = locrb
      locr = locd + n
      locz = locr + n
      locp = locz + n
      locdz = locp + n
      locw  = locdz + n
c
c         check the workspace allocations.
      call dchkw( 'dsdcg', lociw, leniw, locw, lenw, ierr, iter, err )
      if( ierr.ne.0 ) return
c
      iwork(4) = locd
      iwork(9) = lociw
      iwork(10) = locw
c
c         compute the inverse of the diagonal of the matrix.  this
c         will be used as the preconditioner.
      call dsds(n, nelt, ia, ja, a, isym, rwork(locd))
c
c         do the preconditioned conjugate gradient.
      call dcg(n, b, x, nelt, ia, ja, a, isym, dsmv, dsdi,
     $     itol, tol, itmax, iter, err, ierr, iunit, rwork(locr),
     $     rwork(locz), rwork(locp), rwork(locdz), rwork, iwork)
      return
c------------- last line of dsdcg follows -----------------------------
      end
