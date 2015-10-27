*deck domn
      subroutine domn (n, b, x, nelt, ia, ja, a, isym, matvec, msolve,
     +   nsave, itol, tol, itmax, iter, err, ierr, iunit, r, z, p, ap,
     +   emap, dz, csav, rwork, iwork)
c***begin prologue  domn
c***purpose  preconditioned orthomin sparse iterative ax=b solver.
c            routine to solve a general linear system  ax = b  using
c            the preconditioned orthomin method.
c***library   slatec (slap)
c***category  d2a4, d2b4
c***type      double precision (somn-s, domn-d)
c***keywords  iterative precondition, non-symmetric linear system,
c             orthomin, slap, sparse
c***author  greenbaum, anne, (courant institute)
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c *usage:
c     integer  n, nelt, ia(nelt), ja(nelt), isym, nsave, itol, itmax
c     integer  iter, ierr, iunit, iwork(user defined)
c     double precision b(n), x(n), a(nelt), tol, err, r(n), z(n)
c     double precision p(n,0:nsave), ap(n,0:nsave), emap(n,0:nsave)
c     double precision dz(n), csav(nsave), rwork(user defined)
c     external matvec, msolve
c
c     call domn(n, b, x, nelt, ia, ja, a, isym, matvec, msolve,
c    $     nsave, itol, tol, itmax, iter, err, ierr, iunit, r,
c    $     z, p, ap, emap, dz, csav, rwork, iwork)
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
c ia     :in       integer ia(nelt).
c ja     :in       integer ja(nelt).
c a      :in       double precision a(nelt).
c         these arrays contain the matrix data structure for a.
c         it could take any form.  see "description", below, for more
c         details.
c isym   :in       integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all non-zero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the upper
c         or lower triangle of the matrix is stored.
c matvec :ext      external.
c         name of a routine which performs the matrix vector multiply
c         y = a*x given a and x.  the name of the matvec routine must
c         be declared external in the calling program.  the calling
c         sequence to matvec is:
c             call matvec( n, x, y, nelt, ia, ja, a, isym )
c         where n is the number of unknowns, y is the product a*x
c         upon return x is an input vector, nelt is the number of
c         non-zeros in the slap ia, ja, a storage for the matrix a.
c         isym is a flag which, if non-zero, denotest that a is
c         symmetric and only the lower or upper triangle is stored.
c msolve :ext      external.
c         name of a routine which solves a linear system mz = r for
c         z given r with the preconditioning matrix m (m is supplied via
c         rwork and iwork arrays).  the name of the msolve routine must
c         be declared external in the calling program.  the calling
c         sequence to msolve is:
c             call msolve(n, r, z, nelt, ia, ja, a, isym, rwork, iwork)
c         where n is the number of unknowns, r is the right-hand side
c         vector and z is the solution upon return.  nelt, ia, ja, a and
c         isym are defined as above.  rwork is a double precision array
c         that can be used to pass necessary preconditioning information
c         and/or workspace to msolve.  iwork is an integer work array
c         for the same purpose as rwork.
c nsave  :in       integer.
c         number of  direction vectors to save and orthogonalize
c         against.  nsave >= 0.
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
c           ierr = 6 => breakdown of method detected.
c                       (p,ap) < epsilon**2.
c iunit  :in       integer.
c         unit number on which to write the error at each iteration,
c         if this is desired for monitoring convergence.  if unit
c         number is 0, no writing will occur.
c r      :work     double precision r(n).
c z      :work     double precision z(n).
c p      :work     double precision p(n,0:nsave).
c ap     :work     double precision ap(n,0:nsave).
c emap   :work     double precision emap(n,0:nsave).
c dz     :work     double precision dz(n).
c csav   :work     double precision csav(nsave)
c         double precision arrays used for workspace.
c rwork  :work     double precision rwork(user defined).
c         double precision array that can be used for workspace in
c         msolve.
c iwork  :work     integer iwork(user defined).
c         integer array that can be used for workspace in msolve.
c
c *description
c       this routine does  not care  what matrix data   structure is
c       used for  a and m.  it simply   calls  the matvec and msolve
c       routines, with  the arguments as  described above.  the user
c       could write any type of structure and the appropriate matvec
c       and msolve routines.  it is assumed  that a is stored in the
c       ia, ja, a  arrays in some fashion and  that m (or inv(m)) is
c       stored  in  iwork  and  rwork)  in  some fashion.   the slap
c       routines dsdomn and dsluom are examples of this procedure.
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
c       first  in each "column") and are  stored in the  double pre-
c       cision array  a. in  other  words,  for each  column  in the
c       matrix  first put  the diagonal entry in a.  then put in the
c       other non-zero  elements going  down the column  (except the
c       diagonal)  in order.  the ia array  holds the  row index for
c       each non-zero.  the ja array  holds the offsets into the ia,
c       a  arrays  for  the  beginning  of  each  column.  that  is,
c       ia(ja(icol)),a(ja(icol)) are the first elements of the icol-
c       th column in ia and a, and ia(ja(icol+1)-1), a(ja(icol+1)-1)
c       are  the last elements of the icol-th column.   note that we
c       always have ja(n+1)=nelt+1, where n is the number of columns
c       in the matrix  and nelt  is the number  of non-zeros  in the
c       matrix.
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
c***see also  dsdomn, dsluom, isdomn
c***references  1. mark k. seager, a slap for the masses, in
c                  g. f. carey, ed., parallel supercomputing: methods,
c                  algorithms and applications, wiley, 1989, pp.135-155.
c***routines called  d1mach, daxpy, dcopy, ddot, isdomn
c***revision history  (yymmdd)
c   890404  date written
c   890404  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   891004  added new reference.
c   910411  prologue converted to version 4.0 format.  (bab)
c   910502  removed matvec and msolve from routines called list.  (fnf)
c   920407  common block renamed dslblk.  (wrb)
c   920511  added complete declaration section.  (wrb)
c   920929  corrected format of reference.  (fnf)
c   921019  changed 500.0 to 500 to reduce sp/dp differences.  (fnf)
c   921113  corrected c***category line.  (fnf)
c   930326  removed unused variable.  (fnf)
c***end prologue  domn
c     .. scalar arguments ..
      double precision err, tol
      integer ierr, isym, iter, itmax, itol, iunit, n, nelt, nsave
c     .. array arguments ..
      double precision a(nelt), ap(n,0:nsave), b(n), csav(nsave),
     +                 dz(n), emap(n,0:nsave), p(n,0:nsave), r(n),
     +                 rwork(*), x(n), z(n)
      integer ia(nelt), iwork(*), ja(nelt)
c     .. subroutine arguments ..
      external matvec, msolve
c     .. local scalars ..
      double precision ak, akden, aknum, bkl, bnrm, fuzz, solnrm
      integer i, ip, ipo, k, l, lmax
c     .. external functions ..
      double precision d1mach, ddot
      integer isdomn
      external d1mach, ddot, isdomn
c     .. external subroutines ..
      external daxpy, dcopy
c     .. intrinsic functions ..
      intrinsic abs, min, mod
c***first executable statement  domn
c
c         check some of the input data.
c
      iter = 0
      ierr = 0
      if( n.lt.1 ) then
         ierr = 3
         return
      endif
      fuzz = d1mach(3)
      if( tol.lt.500*fuzz ) then
         tol = 500*fuzz
         ierr = 4
      endif
      fuzz = fuzz*fuzz
c
c         calculate initial residual and pseudo-residual, and check
c         stopping criterion.
      call matvec(n, x, r, nelt, ia, ja, a, isym)
      do 10 i = 1, n
         r(i)  = b(i) - r(i)
 10   continue
      call msolve(n, r, z, nelt, ia, ja, a, isym, rwork, iwork)
c
      if( isdomn(n, b, x, nelt, ia, ja, a, isym, msolve, nsave,
     $     itol, tol, itmax, iter, err, ierr, iunit,
     $     r, z, p, ap, emap, dz, csav,
     $     rwork, iwork, ak, bnrm, solnrm) .ne. 0 ) go to 200
      if( ierr.ne.0 ) return
c
c
c         ***** iteration loop *****
c
cvd$r novector
cvd$r noconcur
      do 100 k = 1, itmax
         iter = k
         ip = mod( iter-1, nsave+1 )
c
c         calculate direction vector p, a*p, and (m-inv)*a*p,
c         and save if desired.
         call dcopy(n, z, 1, p(1,ip), 1)
         call matvec(n, p(1,ip), ap(1,ip), nelt, ia, ja, a, isym)
         call msolve(n, ap(1,ip), emap(1,ip), nelt, ia, ja, a, isym,
     $        rwork, iwork)
         if( nsave.eq.0 ) then
            akden = ddot(n, emap, 1, emap, 1)
         else
            if( iter.gt.1 ) then
               lmax = min( nsave, iter-1 )
               do 20 l = 1, lmax
                  ipo = mod(ip+(nsave+1-l),nsave+1)
                  bkl = ddot(n, emap(1,ip), 1, emap(1,ipo), 1)
                  bkl = bkl*csav(l)
                  call daxpy(n, -bkl,    p(1,ipo), 1,    p(1,ip), 1)
                  call daxpy(n, -bkl,   ap(1,ipo), 1,   ap(1,ip), 1)
                  call daxpy(n, -bkl, emap(1,ipo), 1, emap(1,ip), 1)
 20            continue
               if( nsave.gt.1 ) then
                  do 30 l = nsave-1, 1, -1
                     csav(l+1) = csav(l)
 30               continue
               endif
            endif
            akden = ddot(n, emap(1,ip), 1, emap(1,ip), 1)
            if( abs(akden).lt.fuzz ) then
               ierr = 6
               return
            endif
            csav(1) = 1.0d0/akden
c
c         calculate coefficient ak, new iterate x, new residual r, and
c         new pseudo-residual z.
         endif
         aknum = ddot(n, z, 1, emap(1,ip), 1)
         ak = aknum/akden
         call daxpy(n,  ak,    p(1,ip), 1, x, 1)
         call daxpy(n, -ak,   ap(1,ip), 1, r, 1)
         call daxpy(n, -ak, emap(1,ip), 1, z, 1)
c
c         check stopping criterion.
         if( isdomn(n, b, x, nelt, ia, ja, a, isym, msolve, nsave,
     $        itol, tol, itmax, iter, err, ierr, iunit,
     $        r, z, p, ap, emap, dz, csav,
     $        rwork, iwork, ak, bnrm, solnrm) .ne. 0 ) go to 200
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
c------------- last line of domn follows ----------------------------
      end
