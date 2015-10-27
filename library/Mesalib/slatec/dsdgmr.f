*deck dsdgmr
      subroutine dsdgmr (n, b, x, nelt, ia, ja, a, isym, nsave, itol,
     +   tol, itmax, iter, err, ierr, iunit, rwork, lenw, iwork, leniw)
c***begin prologue  dsdgmr
c***purpose  diagonally scaled gmres iterative sparse ax=b solver.
c            this routine uses the generalized minimum residual
c            (gmres) method with diagonal scaling to solve possibly
c            non-symmetric linear systems of the form: ax = b.
c***library   slatec (slap)
c***category  d2a4, d2b4
c***type      double precision (ssdgmr-s, dsdgmr-d)
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
c      integer   n, nelt, ia(nelt), ja(nelt), isym, nsave, itol
c      integer   itmax, iter, ierr, iunit, lenw, iwork(leniw), leniw
c      double precision b(n), x(n), a(nelt), tol, err, rwork(lenw)
c
c      call dsdgmr(n, b, x, nelt, ia, ja, a, isym, nsave,
c     $     itol, tol, itmax, iter, err, ierr, iunit,
c     $     rwork, lenw, iwork, leniw)
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
c         these arrays should hold the matrix a in either the slap
c         triad format or the slap column format.  see "description",
c         below.  if the slap triad format is chosen it is changed
c         internally to the slap column format.
c isym   :in       integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all non-zero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the upper
c         or lower triangle of the matrix is stored.
c nsave  :in       integer.
c         number of direction vectors to save and orthogonalize against.
c         must be greater than 1.
c itol   :in       integer.
c         flag to indicate the type of convergence criterion used.
c         itol=0  means the  iteration stops when the test described
c                 below on  the  residual rl  is satisfied.  this is
c                 the  "natural stopping criteria" for this routine.
c                 other values  of   itol  cause  extra,   otherwise
c                 unnecessary, computation per iteration and     are
c                 therefore  much less  efficient.  see  isdgmr (the
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
c                 there is either left  or no  preconditioning begin
c                 used.
c         itol=11 is    often  useful  for   checking  and comparing
c                 different routines.  for this case, the  user must
c                 supply  the  "exact" solution or  a  very accurate
c                 approximation (one with  an  error much less  than
c                 tol) through a common block,
c                     common /dslblk/ soln( )
c                 if itol=11, iteration stops when the 2-norm of the
c                 difference between the iterative approximation and
c                 the user-supplied solution  divided by the  2-norm
c                 of the  user-supplied solution  is  less than tol.
c                 note that this requires  the  user to  set up  the
c                 "common     /dslblk/ soln(length)"  in the calling
c                 routine.  the routine with this declaration should
c                 be loaded before the stop test so that the correct
c                 length is used by  the loader.  this procedure  is
c                 not standard fortran and may not work correctly on
c                 your   system (although  it  has  worked  on every
c                 system the authors have tried).  if itol is not 11
c                 then this common block is indeed standard fortran.
c tol    :inout    double precision.
c         convergence criterion, as described below.  if tol is set
c         to zero on input, then a default value of 500*(the smallest
c         positive magnitude, machine epsilon) is used.
c itmax  :in       integer.
c         maximum number of iterations.  this routine uses the default
c         of nrmax = itmax/nsave to determine when each restart
c         should occur.  see the description of nrmax and maxl in
c         dgmres for a full and frightfully interesting discussion of
c         this topic.
c iter   :out      integer.
c         number of iterations required to reach convergence, or
c         itmax+1 if convergence criterion could not be achieved in
c         itmax iterations.
c err    :out      double precision.
c         error estimate of error in final approximate solution, as
c         defined by itol.  letting norm() denote the euclidean
c         norm, err is defined as follows...
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
c               ierr = 2 => routine dpigmr failed to reduce the norm
c                           of the current residual on its last call,
c                           and so the iteration has stalled.  in
c                           this case, x equals the last computed
c                           approximation.  the user must either
c                           increase maxl, or choose a different
c                           initial guess.
c               ierr =-1 => insufficient length for rgwk array.
c                           igwk(6) contains the required minimum
c                           length of the rgwk array.
c               ierr =-2 => inconsistent itol and jpre values.
c         for ierr <= 2, rgwk(1) = rhol, which is the norm on the
c         left-hand-side of the relevant stopping test defined
c         below associated with the residual for the current
c         approximation x(l).
c iunit  :in       integer.
c         unit number on which to write the error at each iteration,
c         if this is desired for monitoring convergence.  if unit
c         number is 0, no writing will occur.
c rwork  :work    double precision rwork(lenw).
c         double precision array of size lenw.
c lenw   :in       integer.
c         length of the double precision workspace, rwork.
c         lenw >= 1 + n*(nsave+7) + nsave*(nsave+3).
c         for the recommended values of nsave (10), rwork has size at
c         least 131 + 17*n.
c iwork  :inout    integer iwork(user defined >= 30).
c         used to hold pointers into the rwork array.
c         upon return the following locations of iwork hold information
c         which may be of use to the user:
c         iwork(9)  amount of integer workspace actually used.
c         iwork(10) amount of double precision workspace actually used.
c leniw  :in       integer.
c         length of the integer workspace iwork.  leniw >= 30.
c
c *description:
c       dsdgmr solves a linear system a*x = b rewritten in the form:
c
c        (sb*a*(m-inverse)*(sx-inverse))*(sx*m*x) = sb*b,
c
c       with right preconditioning, or
c
c        (sb*(m-inverse)*a*(sx-inverse))*(sx*x) = sb*(m-inverse)*b,
c
c       with left preconditioning, where a is an n-by-n double precision
c       matrix, x and b are n-vectors, sb and sx are diagonal scaling
c       matrices, and  m   is   the  diagonal  of   a.     it   uses
c       preconditioned   krylov  subpace   methods  based    on  the
c       generalized  minimum residual method (gmres).   this routine
c       is  a  driver routine  which   assumes a  slap matrix   data
c       structure  and   sets  up the  necessary information   to do
c       diagonal preconditioning and  calls  the main gmres  routine
c       dgmres   for  the  solution  of the   linear system.  dgmres
c       optionally   performs   either the   full  orthogonalization
c       version of the gmres algorithm or an  incomplete  variant of
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
c       the user.  if tol equals zero  when dsdgmr is called, then a
c       default  value  of 500*(the   smallest  positive  magnitude,
c       machine epsilon) is used.  if the  scaling arrays sb  and sx
c       are used, then  ideally they  should be chosen  so  that the
c       vectors sx*x(or sx*m*x) and  sb*b have all their  components
c       approximately equal  to  one in  magnitude.  if one wants to
c       use the same scaling in x  and b, then  sb and sx can be the
c       same array in the calling program.
c
c       the following is a list of the other routines and their
c       functions used by gmres:
c       dgmres  contains the matrix structure independent driver
c               routine for gmres.
c       dpigmr  contains the main iteration loop for gmres.
c       dorth   orthogonalizes a new vector against older basis vectors.
c       dheqr   computes a qr decomposition of a hessenberg matrix.
c       dhels   solves a hessenberg least-squares system, using qr
c               factors.
c       rlcalc  computes the scaled residual rl.
c       xlcalc  computes the solution xl.
c       isdgmr  user-replaceable stopping routine.
c
c       the sparse linear algebra package (slap) utilizes two matrix
c       data structures: 1) the  slap triad  format or  2)  the slap
c       column format.  the user can hand this routine either of the
c       of these data structures and slap  will figure out  which on
c       is being used and act accordingly.
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
c***references  1. peter n. brown and a. c. hindmarsh, reduced storage
c                  matrix methods in stiff ode systems, lawrence liver-
c                  more national laboratory report ucrl-95088, rev. 1,
c                  livermore, california, june 1987.
c***routines called  dchkw, dgmres, ds2y, dsdi, dsds, dsmv
c***revision history  (yymmdd)
c   890404  date written
c   890404  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   910411  prologue converted to version 4.0 format.  (bab)
c   920407  common block renamed dslblk.  (wrb)
c   920511  added complete declaration section.  (wrb)
c   920929  corrected format of references.  (fnf)
c***end prologue  dsdgmr
c         the following is for optimized compilation on llnl/ltss crays.
clll. optimize
c     .. parameters ..
      integer locrb, locib
      parameter (locrb=1, locib=11)
c     .. scalar arguments ..
      double precision err, tol
      integer ierr, isym, iter, itmax, itol, iunit, leniw, lenw, n,
     +        nelt, nsave
c     .. array arguments ..
      double precision a(nelt), b(n), rwork(lenw), x(n)
      integer ia(nelt), iwork(leniw), ja(nelt)
c     .. local scalars ..
      integer locdin, locigw, lociw, locrgw, locw, myitol
c     .. external subroutines ..
      external dchkw, dgmres, ds2y, dsdi, dsds, dsmv
c***first executable statement  dsdgmr
c
      ierr = 0
      err  = 0
      if( nsave.le.1 ) then
         ierr = 3
         return
      endif
c
c         change the slap input matrix ia, ja, a to slap-column format.
      call ds2y( n, nelt, ia, ja, a, isym )
c
c         set up the workspace.  we assume maxl=kmp=nsave.
      locigw = locib
      lociw = locigw + 20
c
      locdin = locrb
      locrgw = locdin + n
      locw = locrgw + 1+n*(nsave+6)+nsave*(nsave+3)
c
      iwork(4) = locdin
      iwork(9) = lociw
      iwork(10) = locw
c
c         check the workspace allocations.
      call dchkw( 'dsdgmr', lociw, leniw, locw, lenw, ierr, iter, err )
      if( ierr.ne.0 ) return
c
c         compute the inverse of the diagonal of the matrix.
      call dsds(n, nelt, ia, ja, a, isym, rwork(locdin))
c
c         perform the diagonally scaled generalized minimum
c         residual iteration algorithm.  the following dgmres
c         defaults are used maxl = kmp = nsave, jscal = 0,
c         jpre = -1, nrmax = itmax/nsave
      iwork(locigw  ) = nsave
      iwork(locigw+1) = nsave
      iwork(locigw+2) = 0
      iwork(locigw+3) = -1
      iwork(locigw+4) = itmax/nsave
      myitol = 0
c
      call dgmres( n, b, x, nelt, ia, ja, a, isym, dsmv, dsdi,
     $     myitol, tol, itmax, iter, err, ierr, iunit, rwork, rwork,
     $     rwork(locrgw), lenw-locrgw, iwork(locigw), 20,
     $     rwork, iwork )
c
      if( iter.gt.itmax ) ierr = 2
      return
c------------- last line of dsdgmr follows ----------------------------
      end
