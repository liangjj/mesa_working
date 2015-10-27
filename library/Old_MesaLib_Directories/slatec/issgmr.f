*deck issgmr
      integer function issgmr (n, b, x, xl, nelt, ia, ja, a, isym,
     +   msolve, nmsl, itol, tol, itmax, iter, err, iunit, r, z, dz,
     +   rwork, iwork, rnrm, bnrm, sb, sx, jscal, kmp, lgmr, maxl,
     +   maxlp1, v, q, snormw, prod, r0nrm, hes, jpre)
c***begin prologue  issgmr
c***subsidiary
c***purpose  generalized minimum residual stop test.
c            this routine calculates the stop test for the generalized
c            minimum residual (gmres) iteration scheme.  it returns a
c            non-zero if the error estimate (the type of which is
c            determined by itol) is less than the user specified
c            tolerance tol.
c***library   slatec (slap)
c***category  d2a4, d2b4
c***type      single precision (issgmr-s, isdgmr-d)
c***keywords  gmres, linear system, slap, sparse, stop test
c***author  brown, peter, (llnl), pnbrown@llnl.gov
c           hindmarsh, alan, (llnl), alanh@llnl.gov
c           seager, mark k., (llnl), seager@llnl.gov
c             lawrence livermore national laboratory
c             po box 808, l-60
c             livermore, ca 94550 (510) 423-3141
c***description
c
c *usage:
c      integer n, nelt, ia(nelt), ja(nelt), isym, nmsl, itol
c      integer itmax, iter, iunit, iwork(user defined), jscal
c      integer kmp, lgmr, maxl, maxlp1, jpre
c      real b(n), x(n), xl(maxl), a(nelt), tol, err, r(n), z(n),
c     $     dz(n), rwork(user defined), rnrm, bnrm, sb(n), sx(n),
c     $     v(n,maxlp1), q(2*maxl), snormw, prod, r0nrm,
c     $     hes(maxlp1,maxl)
c      external msolve
c
c      if (issgmr(n, b, x, xl, nelt, ia, ja, a, isym, msolve,
c     $     nmsl, itol, tol, itmax, iter, err, iunit, r, z, dz,
c     $     rwork, iwork, rnrm, bnrm, sb, sx, jscal,
c     $     kmp, lgmr, maxl, maxlp1, v, q, snormw, prod, r0nrm,
c     $     hes, jpre) .ne. 0) then iteration done
c
c *arguments:
c n      :in       integer.
c         order of the matrix.
c b      :in       real b(n).
c         right-hand-side vector.
c x      :in       real x(n).
c         approximate solution vector as of the last restart.
c xl     :out      real xl(n)
c         an array of length n used to hold the approximate
c         solution as of the current iteration.  only computed by
c         this routine when itol=11.
c nelt   :in       integer.
c         number of non-zeros stored in a.
c ia     :in       integer ia(nelt).
c ja     :in       integer ja(nelt).
c a      :in       real a(nelt).
c         these arrays contain the matrix data structure for a.
c         it could take any form.  see "description", in the sgmres,
c         sslugm and ssdgmr routines for more details.
c isym   :in       integer.
c         flag to indicate symmetric storage format.
c         if isym=0, all non-zero entries of the matrix are stored.
c         if isym=1, the matrix is symmetric, and only the upper
c         or lower triangle of the matrix is stored.
c msolve :ext      external.
c         name of a routine which solves a linear system mz = r for  z
c         given r with the preconditioning matrix m (m is supplied via
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
c nmsl   :inout    integer.
c         a counter for the number of calls to msolve.
c itol   :in       integer.
c         flag to indicate the type of convergence criterion used.
c         itol=0  means the  iteration stops when the test described
c                 below on  the  residual rl  is satisfied.  this is
c                 the  "natural stopping criteria" for this routine.
c                 other values  of   itol  cause  extra,   otherwise
c                 unnecessary, computation per iteration and     are
c                 therefore much less efficient.
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
c tol    :in       real.
c         convergence criterion, as described above.
c itmax  :in       integer.
c         maximum number of iterations.
c iter   :in       integer.
c         the iteration for which to check for convergence.
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
c iunit  :in       integer.
c         unit number on which to write the error at each iteration,
c         if this is desired for monitoring convergence.  if unit
c         number is 0, no writing will occur.
c r      :inout    real r(n).
c         work array used in calling routine.  it contains
c         information necessary to compute the residual rl = b-a*xl.
c z      :work     real z(n).
c         workspace used to hold the pseudo-residual m z = r.
c dz     :work     real dz(n).
c         workspace used to hold temporary vector(s).
c rwork  :work     real rwork(user defined).
c         real array that can be used by msolve.
c iwork  :work     integer iwork(user defined).
c         integer array that can be used by msolve.
c rnrm   :in       real.
c         norm of the current residual.  type of norm depends on itol.
c bnrm   :in       real.
c         norm of the right hand side.  type of norm depends on itol.
c sb     :in       real sb(n).
c         scaling vector for b.
c sx     :in       real sx(n).
c         scaling vector for x.
c jscal  :in       integer.
c         flag indicating if scaling arrays sb and sx are being
c         used in the calling routine spigmr.
c         jscal=0 means sb and sx are not used and the
c                 algorithm will perform as if all
c                 sb(i) = 1 and sx(i) = 1.
c         jscal=1 means only sx is used, and the algorithm
c                 performs as if all sb(i) = 1.
c         jscal=2 means only sb is used, and the algorithm
c                 performs as if all sx(i) = 1.
c         jscal=3 means both sb and sx are used.
c kmp    :in       integer
c         the number of previous vectors the new vector vnew
c         must be made orthogonal to.  (kmp .le. maxl)
c lgmr   :in       integer
c         the number of gmres iterations performed on the current call
c         to spigmr (i.e., # iterations since the last restart) and
c         the current order of the upper hessenberg
c         matrix hes.
c maxl   :in       integer
c         the maximum allowable order of the matrix h.
c maxlp1 :in       integer
c         maxpl1 = maxl + 1, used for dynamic dimensioning of hes.
c v      :in       real v(n,maxlp1)
c         the n by (lgmr+1) array containing the lgmr
c         orthogonal vectors v(*,1) to v(*,lgmr).
c q      :in       real q(2*maxl)
c         a real array of length 2*maxl containing the components
c         of the givens rotations used in the qr decomposition
c         of hes.
c snormw :in       real
c         a scalar containing the scaled norm of vnew before it
c         is renormalized in spigmr.
c prod   :in       real
c         the product s1*s2*...*sl = the product of the sines of the
c         givens rotations used in the qr factorization of the
c         hessenberg matrix hes.
c r0nrm  :in       real
c         the scaled norm of initial residual r0.
c hes    :in       real hes(maxlp1,maxl)
c         the upper triangular factor of the qr decomposition
c         of the (lgmr+1) by lgmr upper hessenberg matrix whose
c         entries are the scaled inner-products of a*v(*,i)
c         and v(*,k).
c jpre   :in       integer
c         preconditioner type flag.
c         (see description of igwk(4) in sgmres.)
c
c *description
c       when using the gmres solver,  the preferred value  for itol
c       is 0.  this is due to the fact that when itol=0 the norm of
c       the residual required in the stopping test is  obtained for
c       free, since this value is already  calculated  in the gmres
c       algorithm.   the  variable  rnrm contains the   appropriate
c       norm, which is equal to norm(sb*(rl - a*xl))  when right or
c       no   preconditioning is  being  performed,   and equal   to
c       norm(sb*minv*(rl - a*xl))  when using left preconditioning.
c       here, norm() is the euclidean norm.  nonzero values of itol
c       require  additional work  to  calculate the  actual  scaled
c       residual  or its scaled/preconditioned  form,  and/or   the
c       approximate solution xl.  hence, these values of  itol will
c       not be as efficient as itol=0.
c
c *cautions:
c     this routine will attempt to write to the fortran logical output
c     unit iunit, if iunit .ne. 0.  thus, the user must make sure that
c     this logical unit is attached to a file or terminal before calling
c     this routine with a non-zero value for iunit.  this routine does
c     not check for the validity of a non-zero iunit unit number.
c
c     this routine does not verify that itol has a valid value.
c     the calling routine should make such a test before calling
c     issgmr, as is done in sgmres.
c
c***see also  sgmres
c***routines called  r1mach, scopy, snrm2, srlcal, sscal, sxlcal
c***common blocks    sslblk
c***revision history  (yymmdd)
c   871211  date written
c   881213  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   910411  prologue converted to version 4.0 format.  (bab)
c   910502  corrected conversion errors, etc.  (fnf)
c   910502  removed msolve from routines called list.  (fnf)
c   910506  made subsidiary to sgmres.  (fnf)
c   920407  common block renamed sslblk.  (wrb)
c   920511  added complete declaration section.  (wrb)
c   921113  corrected c***category line.  (fnf)
c***end prologue  issgmr
c     .. scalar arguments ..
      real bnrm, err, prod, r0nrm, rnrm, snormw, tol
      integer isym, iter, itmax, itol, iunit, jpre, jscal, kmp, lgmr,
     +        maxl, maxlp1, n, nelt, nmsl
c     .. array arguments ..
      real a(*), b(*), dz(*), hes(maxlp1, maxl), q(*), r(*), rwork(*),
     +     sb(*), sx(*), v(n,*), x(*), xl(*), z(*)
      integer ia(*), iwork(*), ja(*)
c     .. subroutine arguments ..
      external msolve
c     .. arrays in common ..
      real soln(1)
c     .. local scalars ..
      real dxnrm, fuzz, rat, ratmax, solnrm, tem
      integer i, ielmax
c     .. external functions ..
      real r1mach, snrm2
      external r1mach, snrm2
c     .. external subroutines ..
      external scopy, srlcal, sscal, sxlcal
c     .. intrinsic functions ..
      intrinsic abs, max, sqrt
c     .. common blocks ..
      common /sslblk/ soln
c     .. save statement ..
      save solnrm
c***first executable statement  issgmr
      issgmr = 0
      if ( itol.eq.0 ) then
c
c       use input from spigmr to determine if stop conditions are met.
c
         err = rnrm/bnrm
      endif
      if ( (itol.gt.0) .and. (itol.le.3) ) then
c
c       use srlcal to calculate the scaled residual vector.
c       store answer in r.
c
         if ( lgmr.ne.0 ) call srlcal(n, kmp, lgmr, maxl, v, q, r,
     $                                snormw, prod, r0nrm)
         if ( itol.le.2 ) then
c         err = ||residual||/||righthandside||(2-norms).
            err = snrm2(n, r, 1)/bnrm
c
c         unscale r by r0nrm*prod when kmp < maxl.
c
            if ( (kmp.lt.maxl) .and. (lgmr.ne.0) ) then
               tem = 1.0e0/(r0nrm*prod)
               call sscal(n, tem, r, 1)
            endif
         elseif ( itol.eq.3 ) then
c         err = max |(minv*residual)(i)/x(i)|
c         when jpre .lt. 0, r already contains minv*residual.
            if ( jpre.gt.0 ) then
               call msolve(n, r, dz, nelt, ia, ja, a, isym, rwork,
     $              iwork)
               nmsl = nmsl + 1
            endif
c
c         unscale r by r0nrm*prod when kmp < maxl.
c
            if ( (kmp.lt.maxl) .and. (lgmr.ne.0) ) then
               tem = 1.0e0/(r0nrm*prod)
               call sscal(n, tem, r, 1)
            endif
c
            fuzz = r1mach(1)
            ielmax = 1
            ratmax = abs(dz(1))/max(abs(x(1)),fuzz)
            do 25 i = 2, n
               rat = abs(dz(i))/max(abs(x(i)),fuzz)
               if( rat.gt.ratmax ) then
                  ielmax = i
                  ratmax = rat
               endif
 25         continue
            err = ratmax
            if( ratmax.le.tol ) issgmr = 1
            if( iunit.gt.0 ) write(iunit,1020) iter, ielmax, ratmax
            return
         endif
      endif
      if ( itol.eq.11 ) then
c
c       use sxlcal to calculate the approximate solution xl.
c
         if ( (lgmr.ne.0) .and. (iter.gt.0) ) then
            call sxlcal(n, lgmr, x, xl, xl, hes, maxlp1, q, v, r0nrm,
     $           dz, sx, jscal, jpre, msolve, nmsl, rwork, iwork,
     $           nelt, ia, ja, a, isym)
         elseif ( iter.eq.0 ) then
c         copy x to xl to check if initial guess is good enough.
            call scopy(n, x, 1, xl, 1)
         else
c         return since this is the first call to spigmr on a restart.
            return
         endif
c
         if ((jscal .eq. 0) .or.(jscal .eq. 2)) then
c         err = ||x-truesolution||/||truesolution||(2-norms).
            if ( iter.eq.0 ) solnrm = snrm2(n, soln, 1)
            do 30 i = 1, n
               dz(i) = xl(i) - soln(i)
 30         continue
            err = snrm2(n, dz, 1)/solnrm
         else
            if (iter .eq. 0) then
               solnrm = 0
               do 40 i = 1,n
                  solnrm = solnrm + (sx(i)*soln(i))**2
 40            continue
               solnrm = sqrt(solnrm)
            endif
            dxnrm = 0
            do 50 i = 1,n
               dxnrm = dxnrm + (sx(i)*(xl(i)-soln(i)))**2
 50         continue
            dxnrm = sqrt(dxnrm)
c         err = ||sx*(x-truesolution)||/||sx*truesolution|| (2-norms).
            err = dxnrm/solnrm
         endif
      endif
c
      if( iunit.ne.0 ) then
         if( iter.eq.0 ) then
            write(iunit,1000) n, itol, maxl, kmp
         endif
         write(iunit,1010) iter, rnrm/bnrm, err
      endif
      if ( err.le.tol ) issgmr = 1
c
      return
 1000 format(' generalized minimum residual(',i3,i3,') for ',
     $     'n, itol = ',i5, i5,
     $     /' iter','   natural err est','   error estimate')
 1010 format(1x,i4,1x,e16.7,1x,e16.7)
 1020 format(1x,' iter = ',i5, ' ielmax = ',i5,
     $     ' |r(ielmax)/x(ielmax)| = ',e12.5)
c------------- last line of issgmr follows ----------------------------
      end
