*deck isscg
      integer function isscg (n, b, x, nelt, ia, ja, a, isym, msolve,
     +   itol, tol, itmax, iter, err, ierr, iunit, r, z, p, dz, rwork,
     +   iwork, ak, bk, bnrm, solnrm)
c***begin prologue  isscg
c***subsidiary
c***purpose  preconditioned conjugate gradient stop test.
c            this routine calculates the stop test for the conjugate
c            gradient iteration scheme.  it returns a non-zero if the
c            error estimate (the type of which is determined by itol)
c            is less than the user specified tolerance tol.
c***library   slatec (slap)
c***category  d2b4
c***type      single precision (isscg-s, isdcg-d)
c***keywords  linear system, slap, sparse, stop test
c***author  greenbaum, anne, (courant institute)
c           seager, mark k., (llnl)
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c             seager@llnl.gov
c***description
c
c *usage:
c     integer n, nelt, ia(nelt), ja(nelt), isym, itol, itmax, iter
c     integer ierr, iunit, iwork(user defined)
c     real    b(n), x(n), a(n), tol, err, r(n), z(n)
c     real    p(n), dz(n), rwork(user defined), ak, bk
c     real    bnrm, solnrm
c     external msolve
c
c     if( isscg(n, b, x, nelt, ia, ja, a, isym, msolve, itol, tol,
c    $     itmax, iter, err, ierr, iunit, r, z, p, dz, rwork, iwork,
c    $     ak, bk, bnrm, solnrm) .ne. 0 ) then iteration done
c
c *arguments:
c n      :in       integer.
c         order of the matrix.
c b      :in       real b(n).
c         right-hand side vector.
c x      :in       real x(n).
c         the current approximate solution vector.
c nelt   :in       integer.
c         number of non-zeros stored in a.
c ia     :in       integer ia(nelt).
c ja     :in       integer ja(nelt).
c a      :in       real a(nelt).
c         these arrays should hold the matrix a in either the slap
c         triad format or the slap column format.  see "description"
c         in the scg, ssdcg or ssiccg routines.
c isym   :in       integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all non-zero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the upper
c         or lower triangle of the matrix is stored.
c msolve :ext      external.
c         name of a routine which solves a linear system mz = r for
c         z given r with the preconditioning matrix m (m is supplied via
c         rwork and iwork arrays).  the name of the msolve routine must
c         be declared external in the calling program.  the calling
c         sequence to msolve is:
c             call msolve(n, r, z, nelt, ia, ja, a, isym, rwork, iwork)
c         where n is the number of unknowns, r is the right-hand side
c         vector and z is the solution upon return.  nelt, ia, ja, a and
c         isym are defined as above.  rwork is a real array that can
c         be used to pass necessary preconditioning information and/or
c         workspace to msolve.  iwork is an integer work array for
c         the same purpose as rwork.
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
c tol    :in       real.
c         convergence criterion, as described above.
c itmax  :in       integer.
c         maximum number of iterations.
c iter   :in       integer.
c         current iteration count.  (must be zero on first call.)
c err    :out      real.
c         error estimate of error in the x(n) approximate solution, as
c         defined by itol.
c ierr   :out      integer.
c         error flag.  ierr is set to 3 if itol is not one of the
c         acceptable values, see above.
c iunit  :in       integer.
c         unit number on which to write the error at each iteration,
c         if this is desired for monitoring convergence.  if unit
c         number is 0, no writing will occur.
c r      :in       real r(n).
c         the residual r = b-ax.
c z      :work     real z(n).
c         workspace used to hold the pseudo-residual m z = r.
c p      :in       real p(n).
c         the conjugate direction vector.
c dz     :work     real dz(n).
c         workspace used to hold temporary vector(s).
c rwork  :work     real rwork(user defined).
c         real array that can be used by msolve.
c iwork  :work     integer iwork(user defined).
c         integer array that can be used by msolve.
c ak     :in       real.
c bk     :in       real.
c         current conjugate gradient parameters alpha and beta.
c bnrm   :inout    real.
c         norm of the right hand side.  type of norm depends on itol.
c         calculated only on the first call.
c solnrm :inout    real.
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
c***see also  scg, ssdcg, ssiccg
c***routines called  r1mach, snrm2
c***common blocks    sslblk
c***revision history  (yymmdd)
c   871119  date written
c   881213  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890921  removed tex from comments.  (fnf)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   891003  removed c***refer to line, per mks.
c   910411  prologue converted to version 4.0 format.  (bab)
c   910502  removed msolve from routines called list.  (fnf)
c   910506  made subsidiary to scg.  (fnf)
c   920407  common block renamed sslblk.  (wrb)
c   920511  added complete declaration section.  (wrb)
c   920930  corrected to not print ak,bk when iter=0.  (fnf)
c   921026  changed 1.0e10 to r1mach(2).  (fnf)
c***end prologue  isscg
c     .. scalar arguments ..
      real ak, bk, bnrm, err, solnrm, tol
      integer ierr, isym, iter, itmax, itol, iunit, n, nelt
c     .. array arguments ..
      real a(nelt), b(n), dz(n), p(n), r(n), rwork(*), x(n), z(n)
      integer ia(nelt), iwork(*), ja(nelt)
c     .. subroutine arguments ..
      external msolve
c     .. arrays in common ..
      real soln(1)
c     .. local scalars ..
      integer i
c     .. external functions ..
      real r1mach, snrm2
      external r1mach, snrm2
c     .. common blocks ..
      common /sslblk/ soln
c***first executable statement  isscg
      isscg = 0
c
      if( itol.eq.1 ) then
c         err = ||residual||/||righthandside|| (2-norms).
         if(iter .eq. 0) bnrm = snrm2(n, b, 1)
         err = snrm2(n, r, 1)/bnrm
      else if( itol.eq.2 ) then
c                  -1              -1
c         err = ||m  residual||/||m  righthandside|| (2-norms).
         if(iter .eq. 0) then
            call msolve(n, r, z, nelt, ia, ja, a, isym, rwork, iwork)
            bnrm = snrm2(n, dz, 1)
         endif
         err = snrm2(n, z, 1)/bnrm
      else if( itol.eq.11 ) then
c         err = ||x-truesolution||/||truesolution|| (2-norms).
         if(iter .eq. 0) solnrm = snrm2(n, soln, 1)
         do 10 i = 1, n
            dz(i) = x(i) - soln(i)
 10      continue
         err = snrm2(n, dz, 1)/solnrm
      else
c
c         if we get here itol is not one of the acceptable values.
         err = r1mach(2)
         ierr = 3
      endif
c
      if(iunit .ne. 0) then
         if( iter.eq.0 ) then
            write(iunit,1000) n, itol
            write(iunit,1010) iter, err
         else
            write(iunit,1010) iter, err, ak, bk
         endif
      endif
      if(err .le. tol) isscg = 1
      return
 1000 format(' preconditioned conjugate gradient for ',
     $     'n, itol = ',i5, i5,
     $     /' iter','   error estimate','            alpha',
     $     '             beta')
 1010 format(1x,i4,1x,e16.7,1x,e16.7,1x,e16.7)
c------------- last line of isscg follows ------------------------------
      end
