*deck isdcgs
      integer function isdcgs (n, b, x, nelt, ia, ja, a, isym, matvec,
     +   msolve, itol, tol, itmax, iter, err, ierr, iunit, r, r0, p, q,
     +   u, v1, v2, rwork, iwork, ak, bk, bnrm, solnrm)
c***begin prologue  isdcgs
c***subsidiary
c***purpose  preconditioned biconjugate gradient squared stop test.
c            this routine calculates the stop test for the biconjugate
c            gradient squared iteration scheme.  it returns a non-zero
c            if the error estimate (the type of which is determined by
c            itol) is less than the user specified tolerance tol.
c***library   slatec (slap)
c***category  d2a4, d2b4
c***type      double precision (isscgs-s, isdcgs-d)
c***keywords  iterative precondition, non-symmetric linear system, slap,
c             sparse, stop test
c***author  greenbaum, anne, (courant institute)
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c *usage:
c     integer  n, nelt, ia(nelt), ja(nelt), isym, itol, itmax, iter
c     integer  ierr, iunit, iwork(user defined)
c     double precision b(n), x(n), a(n), tol, err, r(n), r0(n), p(n)
c     double precision q(n), u(n), v1(n), v2(n)
c     double precision rwork(user defined), ak, bk, bnrm, solnrm
c     external matvec, msolve
c
c     if( isdcgs(n, b, x, nelt, ia, ja, a, isym, matvec, msolve, itol,
c    $     tol, itmax, iter, err, ierr, iunit, r, r0, p, q, u, v1,
c    $     v2, rwork, iwork, ak, bk, bnrm, solnrm) .ne. 0 )
c    $     then iteration done
c
c *arguments:
c n      :in       integer
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
c         it could take any form.  see "description" in slap routine
c         dcgs for more details.
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
c         return,  x is an input  vector.  nelt, ia,  ja, a, and  isym
c         define the slap matrix data structure.
c msolve :ext      external.
c         name of a routine which solves a linear system mz = r  for z
c         given r with the preconditioning matrix m (m is supplied via
c         rwork  and iwork arrays).   the name  of  the msolve routine
c         must be declared  external  in the  calling   program.   the
c         calling sequence of msolve is:
c             call msolve(n, r, z, nelt, ia, ja, a, isym, rwork, iwork)
c         where n is the number of unknowns, r is  the right-hand side
c         vector, and z is the solution upon return.  nelt, ia, ja, a,
c         and isym define the slap matrix data structure.
c         rwork is a double precision array that can be used to pass
c         necessary preconditioning information and/or workspace to
c         msolve.
c         iwork is an integer work array for the same purpose as rwork.
c itol   :in       integer.
c         flag to indicate type of convergence criterion.
c         if itol=1, iteration stops when the 2-norm of the residual
c         divided by the 2-norm of the right-hand side is less than tol.
c         this routine must calculate the residual from r = a*x - b.
c         this is unnatural and hence expensive for this type of iter-
c         ative method.  itol=2 is *strongly* recommended.
c         if itol=2, iteration stops when the 2-norm of m-inv times the
c         residual divided by the 2-norm of m-inv times the right hand
c         side is less than tol, where m-inv time a vector is the pre-
c         conditioning step.  this is the *natural* stopping for this
c         iterative method and is *strongly* recommended.
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
c tol    :in       double precision.
c         convergence criterion, as described above.
c itmax  :in       integer.
c         maximum number of iterations.
c iter   :in       integer.
c         current iteration count.  (must be zero on first call.)
c         itmax iterations.
c err    :out      double precision.
c         error estimate of error in final approximate solution, as
c         defined by itol.
c ierr   :out      integer.
c         error flag.  ierr is set to 3 if itol is not one of the
c         acceptable values, see above.
c iunit  :in       integer.
c         unit number on which to write the error at each iteration,
c         if this is desired for monitoring convergence.  if unit
c         number is 0, no writing will occur.
c r      :in       double precision r(n).
c         the residual r = b - ax.
c r0     :work     double precision r0(n).
c p      :dummy    double precision p(n).
c q      :dummy    double precision q(n).
c u      :dummy    double precision u(n).
c v1     :dummy    double precision v1(n).
c         double precision arrays used for workspace.
c v2     :work     double precision v2(n).
c         if itol.eq.1 then v2 is used to hold a * x - b on every call.
c         if itol.eq.2 then v2 is used to hold m-inv * b on the first
c         call.
c         if itol.eq.11 then v2 is used to x - soln.
c rwork  :work     double precision rwork(user defined).
c         double precision array that can be used for workspace in
c         msolve.
c iwork  :work     integer iwork(user defined).
c         integer array that can be used for workspace in msolve.
c ak     :in       double precision.
c         current iterate biconjugate gradient iteration parameter.
c bk     :in       double precision.
c         current iterate biconjugate gradient iteration parameter.
c bnrm   :inout    double precision.
c         norm of the right hand side.  type of norm depends on itol.
c         calculated only on the first call.
c solnrm :inout    double precision.
c         2-norm of the true solution, soln.  only computed and used
c         if itol = 11.
c
c *function return values:
c       0 : error estimate (determined by itol) is *not* less than the
c           specified tolerance, tol.  the iteration must continue.
c       1 : error estimate (determined by itol) is less than the
c           specified tolerance, tol.  the iteration can be considered
c           complete.
c
c *cautions:
c     this routine will attempt to write to the fortran logical output
c     unit iunit, if iunit .ne. 0.  thus, the user must make sure that
c     this logical unit is attached to a file or terminal before calling
c     this routine with a non-zero value for iunit.  this routine does
c     not check for the validity of a non-zero iunit unit number.
c
c***see also  dcgs
c***routines called  d1mach, dnrm2
c***common blocks    dslblk
c***revision history  (yymmdd)
c   890404  date written
c   890404  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   891003  removed c***refer to line, per mks.
c   910411  prologue converted to version 4.0 format.  (bab)
c   910502  removed matvec and msolve from routines called list.  (fnf)
c   910506  made subsidiary to dcgs.  (fnf)
c   920407  common block renamed dslblk.  (wrb)
c   920511  added complete declaration section.  (wrb)
c   920930  corrected to not print ak,bk when iter=0.  (fnf)
c   921026  changed 1.0e10 to d1mach(2) and corrected d to e in
c           output format.  (fnf)
c   921113  corrected c***category line.  (fnf)
c***end prologue  isdcgs
c     .. scalar arguments ..
      double precision ak, bk, bnrm, err, solnrm, tol
      integer ierr, isym, iter, itmax, itol, iunit, n, nelt
c     .. array arguments ..
      double precision a(nelt), b(n), p(n), q(n), r(n), r0(n), rwork(*),
     +                 u(n), v1(n), v2(n), x(n)
      integer ia(nelt), iwork(*), ja(nelt)
c     .. subroutine arguments ..
      external matvec, msolve
c     .. arrays in common ..
      double precision soln(1)
c     .. local scalars ..
      integer i
c     .. external functions ..
      double precision d1mach, dnrm2
      external d1mach, dnrm2
c     .. common blocks ..
      common /dslblk/ soln
c***first executable statement  isdcgs
      isdcgs = 0
c
      if( itol.eq.1 ) then
c         err = ||residual||/||righthandside|| (2-norms).
         if(iter .eq. 0) bnrm = dnrm2(n, b, 1)
         call matvec(n, x, v2, nelt, ia, ja, a, isym )
         do 5 i = 1, n
            v2(i) = v2(i) - b(i)
 5       continue
         err = dnrm2(n, v2, 1)/bnrm
      else if( itol.eq.2 ) then
c                  -1              -1
c         err = ||m  residual||/||m  righthandside|| (2-norms).
         if(iter .eq. 0) then
            call msolve(n, b, v2, nelt, ia, ja, a, isym, rwork, iwork)
            bnrm = dnrm2(n, v2, 1)
         endif
         err = dnrm2(n, r, 1)/bnrm
      else if( itol.eq.11 ) then
c         err = ||x-truesolution||/||truesolution|| (2-norms).
         if(iter .eq. 0) solnrm = dnrm2(n, soln, 1)
         do 10 i = 1, n
            v2(i) = x(i) - soln(i)
 10      continue
         err = dnrm2(n, v2, 1)/solnrm
      else
c
c         if we get here itol is not one of the acceptable values.
         err = d1mach(2)
         ierr = 3
      endif
c
c         print the error and coefficients ak, bk on each step,
c         if desired.
      if(iunit .ne. 0) then
         if( iter.eq.0 ) then
            write(iunit,1000) n, itol
            write(iunit,1010) iter, err
         else
            write(iunit,1010) iter, err, ak, bk
         endif
      endif
      if(err .le. tol) isdcgs = 1
c
      return
 1000 format(' preconditioned biconjugate gradient squared for ',
     $     'n, itol = ',i5, i5,
     $     /' iter','   error estimate','            alpha',
     $     '             beta')
 1010 format(1x,i4,1x,d16.7,1x,d16.7,1x,d16.7)
c------------- last line of isdcgs follows ----------------------------
      end
