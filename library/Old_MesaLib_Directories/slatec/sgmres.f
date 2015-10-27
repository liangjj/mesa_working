*deck sgmres
      subroutine sgmres (n, b, x, nelt, ia, ja, a, isym, matvec, msolve,
     +   itol, tol, itmax, iter, err, ierr, iunit, sb, sx, rgwk, lrgw,
     +   igwk, ligw, rwork, iwork)
c***begin prologue  sgmres
c***purpose  preconditioned gmres iterative sparse ax=b solver.
c            this routine uses the generalized minimum residual
c            (gmres) method with preconditioning to solve
c            non-symmetric linear systems of the form: ax = b.
c***library   slatec (slap)
c***category  d2a4, d2b4
c***type      single precision (sgmres-s, dgmres-d)
c***keywords  generalized minimum residual, iterative precondition,
c             non-symmetric linear system, slap, sparse
c***author  brown, peter, (llnl), pnbrown@llnl.gov
c           hindmarsh, alan, (llnl), alanh@llnl.gov
c           seager, mark k., (llnl), seager@llnl.gov
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c***description
c
c *usage:
c      integer   n, nelt, ia(nelt), ja(nelt), isym, itol, itmax
c      integer   iter, ierr, iunit, lrgw, igwk(ligw), ligw
c      integer   iwork(user defined)
c      real      b(n), x(n), a(nelt), tol, err, sb(n), sx(n)
c      real      rgwk(lrgw), rwork(user defined)
c      external  matvec, msolve
c
c      call sgmres(n, b, x, nelt, ia, ja, a, isym, matvec, msolve,
c     $     itol, tol, itmax, iter, err, ierr, iunit, sb, sx,
c     $     rgwk, lrgw, igwk, ligw, rwork, iwork)
c
c *arguments:
c n      :in       integer.
c         order of the matrix.
c b      :in       real b(n).
c         right-hand side vector.
c x      :inout    real x(n).
c         on input x is your initial guess for the solution vector.
c         on output x is the final approximate solution.
c nelt   :in       integer.
c         number of non-zeros stored in a.
c ia     :in       integer ia(nelt).
c ja     :in       integer ja(nelt).
c a      :in       real a(nelt).
c         these arrays contain the matrix data structure for a.
c         it could take any form.  see "description", below,
c         for more details.
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
c             call matvec(n, x, y, nelt, ia, ja, a, isym)
c         where n is the number of unknowns, y is the product a*x
c         upon return, x is an input vector, and nelt is the number of
c         non-zeros in the slap ia, ja, a storage for the matrix a.
c         isym is a flag which, if non-zero, denotes that a is
c         symmetric and only the lower or upper triangle is stored.
c msolve :ext      external.
c         name of the routine which solves a linear system mz = r for
c         z given r with the preconditioning matrix m (m is supplied via
c         rwork and iwork arrays.  the name of the msolve routine must
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
c         flag to indicate the type of convergence criterion used.
c         itol=0  means the  iteration stops when the test described
c                 below on  the  residual rl  is satisfied.  this is
c                 the  "natural stopping criteria" for this routine.
c                 other values  of   itol  cause  extra,   otherwise
c                 unnecessary, computation per iteration and     are
c                 therefore  much less  efficient.  see  issgmr (the
c                 stop test routine) for more information.
c         itol=1  means   the  iteration stops   when the first test
c                 described below on  the residual rl  is satisfied,
c                 and there  is either right  or  no preconditioning
c                 being used.
c         itol=2  implies     that   the  user    is   using    left
c                 preconditioning, and the second stopping criterion
c                 below is used.
c         itol=3  means the  iteration stops   when  the  third test
c                 described below on minv*residual is satisfied, and
c                 there is either left  or no  preconditioning being
c                 used.
c         itol=11 is    often  useful  for   checking  and comparing
c                 different routines.  for this case, the  user must
c                 supply  the  "exact" solution or  a  very accurate
c                 approximation (one with  an  error much less  than
c                 tol) through a common block,
c                     common /sslblk/ soln( )
c                 if itol=11, iteration stops when the 2-norm of the
c                 difference between the iterative approximation and
c                 the user-supplied solution  divided by the  2-norm
c                 of the  user-supplied solution  is  less than tol.
c                 note that this requires  the  user to  set up  the
c                 "common     /sslblk/ soln(length)"  in the calling
c                 routine.  the routine with this declaration should
c                 be loaded before the stop test so that the correct
c                 length is used by  the loader.  this procedure  is
c                 not standard fortran and may not work correctly on
c                 your   system (although  it  has  worked  on every
c                 system the authors have tried).  if itol is not 11
c                 then this common block is indeed standard fortran.
c tol    :inout    real.
c         convergence criterion, as described below.  if tol is set
c         to zero on input, then a default value of 500*(the smallest
c         positive magnitude, machine epsilon) is used.
c itmax  :dummy    integer.
c         maximum number of iterations in most slap routines.  in
c         this routine this does not make sense.  the maximum number
c         of iterations here is given by itmax = maxl*(nrmax+1).
c         see igwk for definitions of maxl and nrmax.
c iter   :out      integer.
c         number of iterations required to reach convergence, or
c         itmax if convergence criterion could not be achieved in
c         itmax iterations.
c err    :out      real.
c         error estimate of error in final approximate solution, as
c         defined by itol.  letting norm() denote the euclidean
c         norm, err is defined as follows..
c
c         if itol=0, then err = norm(sb*(b-a*x(l)))/norm(sb*b),
c                               for right or no preconditioning, and
c                         err = norm(sb*(m-inverse)*(b-a*x(l)))/
c                                norm(sb*(m-inverse)*b),
c                               for left preconditioning.
c         if itol=1, then err = norm(sb*(b-a*x(l)))/norm(sb*b),
c                               since right or no preconditioning
c                               being used.
c         if itol=2, then err = norm(sb*(m-inverse)*(b-a*x(l)))/
c                                norm(sb*(m-inverse)*b),
c                               since left preconditioning is being
c                               used.
c         if itol=3, then err =  max  |(minv*(b-a*x(l)))(i)/x(i)|
c                               i=1,n
c         if itol=11, then err = norm(sb*(x(l)-soln))/norm(sb*soln).
c ierr   :out      integer.
c         return error flag.
c               ierr = 0 => all went well.
c               ierr = 1 => insufficient storage allocated for
c                           rgwk or igwk.
c               ierr = 2 => routine sgmres failed to reduce the norm
c                           of the current residual on its last call,
c                           and so the iteration has stalled.  in
c                           this case, x equals the last computed
c                           approximation.  the user must either
c                           increase maxl, or choose a different
c                           initial guess.
c               ierr =-1 => insufficient length for rgwk array.
c                           igwk(6) contains the required minimum
c                           length of the rgwk array.
c               ierr =-2 => illegal value of itol, or itol and jpre
c                           values are inconsistent.
c         for ierr <= 2, rgwk(1) = rhol, which is the norm on the
c         left-hand-side of the relevant stopping test defined
c         below associated with the residual for the current
c         approximation x(l).
c iunit  :in       integer.
c         unit number on which to write the error at each iteration,
c         if this is desired for monitoring convergence.  if unit
c         number is 0, no writing will occur.
c sb     :in       real sb(n).
c         array of length n containing scale factors for the right
c         hand side vector b.  if jscal.eq.0 (see below), sb need
c         not be supplied.
c sx     :in       real sx(n).
c         array of length n containing scale factors for the solution
c         vector x.  if jscal.eq.0 (see below), sx need not be
c         supplied.  sb and sx can be the same array in the calling
c         program if desired.
c rgwk   :inout    real rgwk(lrgw).
c         real array used for workspace by sgmres.
c         on return, rgwk(1) = rhol.  see ierr for definition of rhol.
c lrgw   :in       integer.
c         length of the real workspace, rgwk.
c         lrgw >= 1 + n*(maxl+6) + maxl*(maxl+3).
c         see below for definition of maxl.
c         for the default values, rgwk has size at least 131 + 16*n.
c igwk   :inout    integer igwk(ligw).
c         the following igwk parameters should be set by the user
c         before calling this routine.
c         igwk(1) = maxl.  maximum dimension of krylov subspace in
c            which x - x0 is to be found (where, x0 is the initial
c            guess).  the default value of maxl is 10.
c         igwk(2) = kmp.  maximum number of previous krylov basis
c            vectors to which each new basis vector is made orthogonal.
c            the default value of kmp is maxl.
c         igwk(3) = jscal.  flag indicating whether the scaling
c            arrays sb and sx are to be used.
c            jscal = 0 => sb and sx are not used and the algorithm
c               will perform as if all sb(i) = 1 and sx(i) = 1.
c            jscal = 1 =>  only sx is used, and the algorithm
c               performs as if all sb(i) = 1.
c            jscal = 2 =>  only sb is used, and the algorithm
c               performs as if all sx(i) = 1.
c            jscal = 3 =>  both sb and sx are used.
c         igwk(4) = jpre.  flag indicating whether preconditioning
c            is being used.
c            jpre = 0  =>  there is no preconditioning.
c            jpre > 0  =>  there is preconditioning on the right
c               only, and the solver will call routine msolve.
c            jpre < 0  =>  there is preconditioning on the left
c               only, and the solver will call routine msolve.
c         igwk(5) = nrmax.  maximum number of restarts of the
c            krylov iteration.  the default value of nrmax = 10.
c            if iwork(5) = -1,  then no restarts are performed (in
c            this case, nrmax is set to zero internally).
c         the following iwork parameters are diagnostic information
c         made available to the user after this routine completes.
c         igwk(6) = mlwk.  required minimum length of rgwk array.
c         igwk(7) = nms.  the total number of calls to msolve.
c ligw   :in       integer.
c         length of the integer workspace, igwk.  ligw >= 20.
c rwork  :work     real rwork(user defined).
c         real array that can be used for workspace in msolve.
c iwork  :work     integer iwork(user defined).
c         integer array that can be used for workspace in msolve.
c
c *description:
c       sgmres solves a linear system a*x = b rewritten in the form:
c
c        (sb*a*(m-inverse)*(sx-inverse))*(sx*m*x) = sb*b,
c
c       with right preconditioning, or
c
c        (sb*(m-inverse)*a*(sx-inverse))*(sx*x) = sb*(m-inverse)*b,
c
c       with left preconditioning, where a is an n-by-n real matrix,
c       x  and  b are n-vectors,   sb and sx   are  diagonal scaling
c       matrices,   and m is  a preconditioning    matrix.   it uses
c       preconditioned  krylov   subpace  methods  based     on  the
c       generalized minimum residual  method (gmres).   this routine
c       optionally performs  either  the  full     orthogonalization
c       version of the  gmres  algorithm or an incomplete variant of
c       it.  both versions use restarting of the linear iteration by
c       default, although the user can disable this feature.
c
c       the gmres  algorithm generates a sequence  of approximations
c       x(l) to the  true solution of the above  linear system.  the
c       convergence criteria for stopping the  iteration is based on
c       the size  of the  scaled norm of  the residual  r(l)  =  b -
c       a*x(l).  the actual stopping test is either:
c
c               norm(sb*(b-a*x(l))) .le. tol*norm(sb*b),
c
c       for right preconditioning, or
c
c               norm(sb*(m-inverse)*(b-a*x(l))) .le.
c                       tol*norm(sb*(m-inverse)*b),
c
c       for left preconditioning, where norm() denotes the euclidean
c       norm, and tol is  a positive scalar less  than one  input by
c       the user.  if tol equals zero  when sgmres is called, then a
c       default  value  of 500*(the   smallest  positive  magnitude,
c       machine epsilon) is used.  if the  scaling arrays sb  and sx
c       are used, then  ideally they  should be chosen  so  that the
c       vectors sx*x(or sx*m*x) and  sb*b have all their  components
c       approximately equal  to  one in  magnitude.  if one wants to
c       use the same scaling in x  and b, then  sb and sx can be the
c       same array in the calling program.
c
c       the following is a list of the other routines and their
c       functions used by sgmres:
c       spigmr  contains the main iteration loop for gmres.
c       sorth   orthogonalizes a new vector against older basis vectors.
c       sheqr   computes a qr decomposition of a hessenberg matrix.
c       shels   solves a hessenberg least-squares system, using qr
c               factors.
c       srlcal  computes the scaled residual rl.
c       sxlcal  computes the solution xl.
c       issgmr  user-replaceable stopping routine.
c
c       this routine does  not care  what matrix data   structure is
c       used for  a and m.  it simply   calls  the matvec and msolve
c       routines, with  the arguments as  described above.  the user
c       could write any type of structure and the appropriate matvec
c       and msolve routines.  it is assumed  that a is stored in the
c       ia, ja, a  arrays in some fashion and  that m (or inv(m)) is
c       stored  in  iwork  and  rwork   in  some fashion.   the slap
c       routines ssdcg and ssiccg are examples of this procedure.
c
c       two  examples  of  matrix  data structures  are the: 1) slap
c       triad  format and 2) slap column format.
c
c       =================== s l a p triad format ===================
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
c *cautions:
c     this routine will attempt to write to the fortran logical output
c     unit iunit, if iunit .ne. 0.  thus, the user must make sure that
c     this logical unit is attached to a file or terminal before calling
c     this routine with a non-zero value for iunit.  this routine does
c     not check for the validity of a non-zero iunit unit number.
c
c***references  1. peter n. brown and a. c. hindmarsh, reduced storage
c                  matrix methods in stiff ode systems, lawrence liver-
c                  more national laboratory report ucrl-95088, rev. 1,
c                  livermore, california, june 1987.
c               2. mark k. seager, a slap for the masses, in
c                  g. f. carey, ed., parallel supercomputing: methods,
c                  algorithms and applications, wiley, 1989, pp.135-155.
c***routines called  r1mach, scopy, snrm2, spigmr
c***revision history  (yymmdd)
c   871001  date written
c   881213  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   891004  added new reference.
c   910411  prologue converted to version 4.0 format.  (bab)
c   910506  corrected errors in c***routines called list.  (fnf)
c   920407  common block renamed sslblk.  (wrb)
c   920511  added complete declaration section.  (wrb)
c   920929  corrected format of references.  (fnf)
c   921019  changed 500.0 to 500 to reduce sp/dp differences.  (fnf)
c   921026  added check for valid value of itol.  (fnf)
c***end prologue  sgmres
c         the following is for optimized compilation on llnl/ltss crays.
clll. optimize
c     .. scalar arguments ..
      real err, tol
      integer ierr, isym, iter, itmax, itol, iunit, ligw, lrgw, n, nelt
c     .. array arguments ..
      real a(nelt), b(n), rgwk(lrgw), rwork(*), sb(n), sx(n), x(n)
      integer ia(nelt), igwk(ligw), iwork(*), ja(nelt)
c     .. subroutine arguments ..
      external matvec, msolve
c     .. local scalars ..
      real bnrm, rhol, sum
      integer i, iflag, jpre, jscal, kmp, ldl, lgmr, lhes, lq, lr, lv,
     +        lw, lxl, lz, lzm1, maxl, maxlp1, nms, nmsl, nrmax, nrsts
c     .. external functions ..
      real r1mach, snrm2
      external r1mach, snrm2
c     .. external subroutines ..
      external scopy, spigmr
c     .. intrinsic functions ..
      intrinsic sqrt
c***first executable statement  sgmres
      ierr = 0
c   ------------------------------------------------------------------
c         load method parameters with user values or defaults.
c   ------------------------------------------------------------------
      maxl = igwk(1)
      if (maxl .eq. 0) maxl = 10
      if (maxl .gt. n) maxl = n
      kmp = igwk(2)
      if (kmp .eq. 0) kmp = maxl
      if (kmp .gt. maxl) kmp = maxl
      jscal = igwk(3)
      jpre = igwk(4)
c         check for valid value of itol.
      if( (itol.lt.0) .or. ((itol.gt.3).and.(itol.ne.11)) ) goto 650
c         check for consistent values of itol and jpre.
      if( itol.eq.1 .and. jpre.lt.0 ) goto 650
      if( itol.eq.2 .and. jpre.ge.0 ) goto 650
      nrmax = igwk(5)
      if( nrmax.eq.0 ) nrmax = 10
c         if nrmax .eq. -1, then set nrmax = 0 to turn off restarting.
      if( nrmax.eq.-1 ) nrmax = 0
c         if input value of tol is zero, set it to its default value.
      if( tol.eq.0.0e0 ) tol = 500*r1mach(3)
c
c         initialize counters.
      iter = 0
      nms = 0
      nrsts = 0
c   ------------------------------------------------------------------
c         form work array segment pointers.
c   ------------------------------------------------------------------
      maxlp1 = maxl + 1
      lv = 1
      lr = lv + n*maxlp1
      lhes = lr + n + 1
      lq = lhes + maxl*maxlp1
      ldl = lq + 2*maxl
      lw = ldl + n
      lxl = lw + n
      lz = lxl + n
c
c         load igwk(6) with required minimum length of the rgwk array.
      igwk(6) = lz + n - 1
      if( lz+n-1.gt.lrgw ) goto 640
c   ------------------------------------------------------------------
c         calculate scaled-preconditioned norm of rhs vector b.
c   ------------------------------------------------------------------
      if (jpre .lt. 0) then
         call msolve(n, b, rgwk(lr), nelt, ia, ja, a, isym,
     $        rwork, iwork)
         nms = nms + 1
      else
         call scopy(n, b, 1, rgwk(lr), 1)
      endif
      if( jscal.eq.2 .or. jscal.eq.3 ) then
         sum = 0
         do 10 i = 1,n
            sum = sum + (rgwk(lr-1+i)*sb(i))**2
 10      continue
         bnrm = sqrt(sum)
      else
         bnrm = snrm2(n,rgwk(lr),1)
      endif
c   ------------------------------------------------------------------
c         calculate initial residual.
c   ------------------------------------------------------------------
      call matvec(n, x, rgwk(lr), nelt, ia, ja, a, isym)
      do 50 i = 1,n
         rgwk(lr-1+i) = b(i) - rgwk(lr-1+i)
 50   continue
c   ------------------------------------------------------------------
c         if performing restarting, then load the residual into the
c         correct location in the rgwk array.
c   ------------------------------------------------------------------
 100  continue
      if( nrsts.gt.nrmax ) goto 610
      if( nrsts.gt.0 ) then
c         copy the current residual to a different location in the rgwk
c         array.
         call scopy(n, rgwk(ldl), 1, rgwk(lr), 1)
      endif
c   ------------------------------------------------------------------
c         use the spigmr algorithm to solve the linear system a*z = r.
c   ------------------------------------------------------------------
      call spigmr(n, rgwk(lr), sb, sx, jscal, maxl, maxlp1, kmp,
     $       nrsts, jpre, matvec, msolve, nmsl, rgwk(lz), rgwk(lv),
     $       rgwk(lhes), rgwk(lq), lgmr, rwork, iwork, rgwk(lw),
     $       rgwk(ldl), rhol, nrmax, b, bnrm, x, rgwk(lxl), itol,
     $       tol, nelt, ia, ja, a, isym, iunit, iflag, err)
      iter = iter + lgmr
      nms = nms + nmsl
c
c         increment x by the current approximate solution z of a*z = r.
c
      lzm1 = lz - 1
      do 110 i = 1,n
         x(i) = x(i) + rgwk(lzm1+i)
 110  continue
      if( iflag.eq.0 ) goto 600
      if( iflag.eq.1 ) then
         nrsts = nrsts + 1
         goto 100
      endif
      if( iflag.eq.2 ) goto 620
c   ------------------------------------------------------------------
c         all returns are made through this section.
c   ------------------------------------------------------------------
c         the iteration has converged.
c
 600  continue
      igwk(7) = nms
      rgwk(1) = rhol
      ierr = 0
      return
c
c         max number((nrmax+1)*maxl) of linear iterations performed.
 610  continue
      igwk(7) = nms
      rgwk(1) = rhol
      ierr = 1
      return
c
c         gmres failed to reduce last residual in maxl iterations.
c         the iteration has stalled.
 620  continue
      igwk(7) = nms
      rgwk(1) = rhol
      ierr = 2
      return
c         error return.  insufficient length for rgwk array.
 640  continue
      err = tol
      ierr = -1
      return
c         error return.  inconsistent itol and jpre values.
 650  continue
      err = tol
      ierr = -2
      return
c------------- last line of sgmres follows ----------------------------
      end
