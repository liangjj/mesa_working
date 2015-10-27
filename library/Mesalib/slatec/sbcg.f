*deck sbcg
      subroutine sbcg (n, b, x, nelt, ia, ja, a, isym, matvec, mttvec,
     +   msolve, mtsolv, itol, tol, itmax, iter, err, ierr, iunit, r, z,
     +   p, rr, zz, pp, dz, rwork, iwork)
c***begin prologue  sbcg
c***purpose  preconditioned biconjugate gradient sparse ax = b solver.
c            routine to solve a non-symmetric linear system  ax = b
c            using the preconditioned biconjugate gradient method.
c***library   slatec (slap)
c***category  d2a4, d2b4
c***type      single precision (sbcg-s, dbcg-d)
c***keywords  biconjugate gradient, iterative precondition,
c             non-symmetric linear system, slap, sparse
c***author  greenbaum, anne, (courant institute)
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c *usage:
c      integer n, nelt, ia(nelt), ja(nelt), isym, itol, itmax
c      integer iter, ierr, iunit, iwork(user defined)
c      real    b(n), x(n), a(nelt), tol, err, r(n), z(n), p(n)
c      real    rr(n), zz(n), pp(n), dz(n)
c      real    rwork(user defined)
c      external matvec, mttvec, msolve, mtsolv
c
c      call sbcg(n, b, x, nelt, ia, ja, a, isym, matvec, mttvec,
c     $     msolve, mtsolv, itol, tol, itmax, iter, err, ierr, iunit,
c     $     r, z, p, rr, zz, pp, dz, rwork, iwork)
c
c *arguments:
c n      :in       integer
c         order of the matrix.
c b      :in       real b(n).
c         right-hand side vector.
c x      :inout    real x(n).
c         on input x is your initial guess for solution vector.
c         on output x is the final approximate solution.
c nelt   :in       integer.
c         number of non-zeros stored in a.
c ia     :in       integer ia(nelt).
c ja     :in       integer ja(nelt).
c a      :in       real a(nelt).
c         these arrays contain the matrix data structure for a.
c         it could take any form.  see "description", below, for more
c         details.
c isym   :in       integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all non-zero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the upper
c         or lower triangle of the matrix is stored.
c matvec :ext      external.
c         name of a routine which  performs the matrix vector multiply
c         operation  y = a*x  given a and x.  the  name of  the matvec
c         routine must  be declared external  in the  calling program.
c         the calling sequence of matvec is:
c             call matvec( n, x, y, nelt, ia, ja, a, isym )
c         where n is the number of unknowns, y is the product a*x upon
c         return,  x is an input  vector.  nelt, ia,  ja,  a and  isym
c         define the slap matrix data structure: see description,below.
c mttvec :ext      external.
c         name of a routine which performs the matrix transpose vector
c         multiply y = a'*x given a and x (where ' denotes transpose).
c         the name of the mttvec routine must be declared external  in
c         the calling program.  the calling sequence to mttvec is  the
c         same as that for mttvec, viz.:
c             call mttvec( n, x, y, nelt, ia, ja, a, isym )
c         where n  is the number  of unknowns, y is the   product a'*x
c         upon return, x is an input vector.  nelt, ia, ja, a and isym
c         define the slap matrix data structure: see description,below.
c msolve :ext      external.
c         name of a routine which solves a linear system mz = r  for z
c         given r with the preconditioning matrix m (m is supplied via
c         rwork  and iwork arrays).   the name  of  the msolve routine
c         must be declared  external  in the  calling   program.   the
c         calling sequence of msolve is:
c             call msolve(n, r, z, nelt, ia, ja, a, isym, rwork, iwork)
c         where n is the number of unknowns, r is  the right-hand side
c         vector, and z is the solution upon return.  nelt,  ia, ja, a
c         and  isym define the slap  matrix  data structure: see
c         description, below.  rwork is a  real array that can be used
c         to  pass   necessary  preconditioning     information and/or
c         workspace to msolve.  iwork is an integer work array for the
c         same purpose as rwork.
c mtsolv :ext      external.
c         name of a routine which solves a linear system m'zz = rr for
c         zz given rr with the preconditioning matrix m (m is supplied
c         via rwork and iwork arrays).  the name of the mtsolv routine
c         must be declared external in the calling program.  the call-
c         ing sequence to mtsolv is:
c            call mtsolv(n, rr, zz, nelt, ia, ja, a, isym, rwork, iwork)
c         where n is the number of unknowns, rr is the right-hand side
c         vector, and zz is the solution upon return.  nelt, ia, ja, a
c         and  isym define the slap  matrix  data structure: see
c         description, below.  rwork is a  real array that can be used
c         to  pass   necessary  preconditioning     information and/or
c         workspace to mtsolv.  iwork is an integer work array for the
c         same purpose as rwork.
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
c r      :work     real r(n).
c z      :work     real z(n).
c p      :work     real p(n).
c rr     :work     real rr(n).
c zz     :work     real zz(n).
c pp     :work     real pp(n).
c dz     :work     real dz(n).
c         real arrays used for workspace.
c rwork  :work     real rwork(user defined).
c         real array that can be used for workspace in msolve
c         and mtsolv.
c iwork  :work     integer iwork(user defined).
c         integer array that can be used for workspace in msolve
c         and mtsolv.
c
c *description
c      this routine does not care what matrix data structure is used
c       for a and m.  it simply calls matvec, mttvec, msolve, mtsolv
c       routines, with arguments as above.  the user could write any
c       type of structure, and  appropriate  matvec, msolve, mttvec,
c       and mtsolv routines.  it  is assumed that a is stored in the
c       ia, ja, a  arrays in some fashion and  that m (or inv(m)) is
c       stored  in  iwork  and  rwork   in  some fashion.   the slap
c       routines ssdbcg and sslubc are examples of this procedure.
c
c       two  examples  of  matrix  data structures  are the: 1) slap
c       triad  format and 2) slap column format.
c
c       =================== s l a p triad format ===================
c       in  this   format only the  non-zeros are  stored.  they may
c       appear  in *any* order.   the user  supplies three arrays of
c       length nelt, where  nelt  is the number  of non-zeros in the
c       matrix:  (ia(nelt), ja(nelt),  a(nelt)).  for each  non-zero
c       the  user puts   the row  and  column index   of that matrix
c       element in the ia and ja arrays.  the  value of the non-zero
c       matrix  element is  placed in  the corresponding location of
c       the a  array.  this is  an extremely easy data  structure to
c       generate.  on  the other hand it  is  not too  efficient  on
c       vector  computers   for the  iterative  solution  of  linear
c       systems.  hence, slap  changes this input  data structure to
c       the slap   column  format for the  iteration (but   does not
c       change it back).
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
c       in  this format   the non-zeros are    stored counting  down
c       columns (except  for the diagonal  entry, which must  appear
c       first in each "column") and are  stored in the real array a.
c       in other words,  for  each column    in the matrix   put the
c       diagonal  entry  in a.   then   put  in the  other  non-zero
c       elements going   down the  column (except  the  diagonal) in
c       order.  the ia array holds the row index  for each non-zero.
c       the ja array holds the offsets into the ia, a arrays for the
c       beginning   of   each  column.      that is,   ia(ja(icol)),
c       a(ja(icol)) points to the beginning of the icol-th column in
c       ia and  a.  ia(ja(icol+1)-1), a(ja(icol+1)-1)  points to the
c       end of the icol-th column.  note that we always have ja(n+1)
c       = nelt+1, where n is the number of columns in the matrix and
c       nelt is the number of non-zeros in the matrix.
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
c *cautions:
c     this routine will attempt to write to the fortran logical output
c     unit iunit, if iunit .ne. 0.  thus, the user must make sure that
c     this logical unit is attached to a file or terminal before calling
c     this routine with a non-zero value for iunit.  this routine does
c     not check for the validity of a non-zero iunit unit number.
c
c***see also  ssdbcg, sslubc
c***references  1. mark k. seager, a slap for the masses, in
c                  g. f. carey, ed., parallel supercomputing: methods,
c                  algorithms and applications, wiley, 1989, pp.135-155.
c***routines called  issbcg, r1mach, saxpy, scopy, sdot
c***revision history  (yymmdd)
c   871119  date written
c   881213  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890921  removed tex from comments.  (fnf)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   891004  added new reference.
c   910411  prologue converted to version 4.0 format.  (bab)
c   910502  removed matvec, mttvec, msolve, mtsolv from routines
c           called list.  (fnf)
c   920407  common block renamed sslblk.  (wrb)
c   920511  added complete declaration section.  (wrb)
c   920929  corrected format of reference.  (fnf)
c   921019  changed 500.0 to 500 to reduce sp/dp differences.  (fnf)
c   921113  corrected c***category line.  (fnf)
c***end prologue  sbcg
c     .. scalar arguments ..
      real err, tol
      integer ierr, isym, iter, itmax, itol, iunit, n, nelt
c     .. array arguments ..
      real a(nelt), b(n), dz(n), p(n), pp(n), r(n), rr(n), rwork(*),
     +     x(n), z(n), zz(n)
      integer ia(nelt), iwork(*), ja(nelt)
c     .. subroutine arguments ..
      external matvec, msolve, mtsolv, mttvec
c     .. local scalars ..
      real ak, akden, bk, bkden, bknum, bnrm, fuzz, solnrm, tolmin
      integer i, k
c     .. external functions ..
      real r1mach, sdot
      integer issbcg
      external r1mach, sdot, issbcg
c     .. external subroutines ..
      external saxpy, scopy
c     .. intrinsic functions ..
      intrinsic abs
c***first executable statement  sbcg
c
c         check some of the input data.
c
      iter = 0
      ierr = 0
      if( n.lt.1 ) then
         ierr = 3
         return
      endif
      fuzz = r1mach(3)
      tolmin = 500*fuzz
      fuzz = fuzz*fuzz
      if( tol.lt.tolmin ) then
         tol = tolmin
         ierr = 4
      endif
c
c         calculate initial residual and pseudo-residual, and check
c         stopping criterion.
      call matvec(n, x, r, nelt, ia, ja, a, isym)
      do 10 i = 1, n
         r(i)  = b(i) - r(i)
         rr(i) = r(i)
 10   continue
      call msolve(n, r, z, nelt, ia, ja, a, isym, rwork, iwork)
      call mtsolv(n, rr, zz, nelt, ia, ja, a, isym, rwork, iwork)
c
      if( issbcg(n, b, x, nelt, ia, ja, a, isym, msolve, itol, tol,
     $     itmax, iter, err, ierr, iunit, r, z, p, rr, zz, pp,
     $     dz, rwork, iwork, ak, bk, bnrm, solnrm) .ne. 0 )
     $     go to 200
      if( ierr.ne.0 ) return
c
c         ***** iteration loop *****
c
      do 100 k=1,itmax
         iter = k
c
c         calculate coefficient bk and direction vectors p and pp.
         bknum = sdot(n, z, 1, rr, 1)
         if( abs(bknum).le.fuzz ) then
            ierr = 6
            return
         endif
         if(iter .eq. 1) then
            call scopy(n, z, 1, p, 1)
            call scopy(n, zz, 1, pp, 1)
         else
            bk = bknum/bkden
            do 20 i = 1, n
               p(i) = z(i) + bk*p(i)
               pp(i) = zz(i) + bk*pp(i)
 20         continue
         endif
         bkden = bknum
c
c         calculate coefficient ak, new iterate x, new residuals r and
c         rr, and new pseudo-residuals z and zz.
         call matvec(n, p, z, nelt, ia, ja, a, isym)
         akden = sdot(n, pp, 1, z, 1)
         ak = bknum/akden
         if( abs(akden).le.fuzz ) then
            ierr = 6
            return
         endif
         call saxpy(n, ak, p, 1, x, 1)
         call saxpy(n, -ak, z, 1, r, 1)
         call mttvec(n, pp, zz, nelt, ia, ja, a, isym)
         call saxpy(n, -ak, zz, 1, rr, 1)
         call msolve(n, r, z, nelt, ia, ja, a, isym, rwork, iwork)
         call mtsolv(n, rr, zz, nelt, ia, ja, a, isym, rwork, iwork)
c
c         check stopping criterion.
         if( issbcg(n, b, x, nelt, ia, ja, a, isym, msolve, itol, tol,
     $        itmax, iter, err, ierr, iunit, r, z, p, rr, zz,
     $        pp, dz, rwork, iwork, ak, bk, bnrm, solnrm) .ne. 0 )
     $        go to 200
c
 100  continue
c
c         *****   end of loop  *****
c
c         stopping criterion not satisfied.
      iter = itmax + 1
      ierr = 2
c
 200  return
c------------- last line of sbcg follows ----------------------------
      end
