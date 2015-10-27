*deck hstcsp
      subroutine hstcsp (intl, a, b, m, mbdcnd, bda, bdb, c, d, n,
     +   nbdcnd, bdc, bdd, elmbda, f, idimf, pertrb, ierror, w)
c***begin prologue  hstcsp
c***purpose  solve the standard five-point finite difference
c            approximation on a staggered grid to the modified helmholtz
c            equation in spherical coordinates assuming axisymmetry
c            (no dependence on longitude).
c***library   slatec (fishpack)
c***category  i2b1a1a
c***type      single precision (hstcsp-s)
c***keywords  elliptic, fishpack, helmholtz, pde, spherical
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c     hstcsp solves the standard five-point finite difference
c     approximation on a staggered grid to the modified helmholtz
c     equation spherical coordinates assuming axisymmetry (no dependence
c     on longitude).
c
c                  (1/r**2)(d/dr)(r**2(du/dr)) +
c
c       1/(r**2*sin(theta))(d/dtheta)(sin(theta)(du/dtheta)) +
c
c            (lambda/(r*sin(theta))**2)u  =  f(theta,r)
c
c     where theta is colatitude and r is the radial coordinate.
c     this two-dimensional modified helmholtz equation results from
c     the fourier transform of the three-dimensional poisson equation.
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c    * * * * * * * *    parameter description     * * * * * * * * * *
c
c
c            * * * * * *   on input    * * * * * *
c
c   intl
c     = 0  on initial entry to hstcsp or if any of the arguments
c          c, d, n, or nbdcnd are changed from a previous call.
c
c     = 1  if c, d, n, and nbdcnd are all unchanged from previous
c          call to hstcsp.
c
c     note:  a call with intl = 0 takes approximately 1.5 times as much
c            time as a call with intl = 1.  once a call with intl = 0
c            has been made then subsequent solutions corresponding to
c            different f, bda, bdb, bdc, and bdd can be obtained
c            faster with intl = 1 since initialization is not repeated.
c
c   a,b
c     the range of theta (colatitude), i.e. a .le. theta .le. b.  a
c     must be less than b and a must be non-negative.  a and b are in
c     radians.  a = 0 corresponds to the north pole and b = pi
c     corresponds to the south pole.
c
c                  * * *  important  * * *
c
c     if b is equal to pi, then b must be computed using the statement
c
c     b = pimach(dum)
c
c     this insures that b in the user's program is equal to pi in this
c     program which permits several tests of the input parameters that
c     otherwise would not be possible.
c
c                  * * * * * * * * * * * *
c
c   m
c     the number of grid points in the interval (a,b).  the grid points
c     in the theta-direction are given by theta(i) = a + (i-0.5)dtheta
c     for i=1,2,...,m where dtheta =(b-a)/m.  m must be greater than 4.
c
c   mbdcnd
c     indicates the type of boundary conditions at theta = a and
c     theta = b.
c
c     = 1  if the solution is specified at theta = a and theta = b.
c          (see notes 1, 2 below)
c
c     = 2  if the solution is specified at theta = a and the derivative
c          of the solution with respect to theta is specified at
c          theta = b (see notes 1, 2 below).
c
c     = 3  if the derivative of the solution with respect to theta is
c          specified at theta = a (see notes 1, 2 below) and theta = b.
c
c     = 4  if the derivative of the solution with respect to theta is
c          specified at theta = a (see notes 1, 2 below) and the
c          solution is specified at theta = b.
c
c     = 5  if the solution is unspecified at theta = a = 0 and the
c          solution is specified at theta = b. (see note 2 below)
c
c     = 6  if the solution is unspecified at theta = a = 0 and the
c          derivative of the solution with respect to theta is
c          specified at theta = b (see note 2 below).
c
c     = 7  if the solution is specified at theta = a and the
c          solution is unspecified at theta = b = pi.
c
c     = 8  if the derivative of the solution with respect to
c          theta is specified at theta = a (see note 1 below)
c          and the solution is unspecified at theta = b = pi.
c
c     = 9  if the solution is unspecified at theta = a = 0 and
c          theta = b = pi.
c
c     notes:  1.  if a = 0, do not use mbdcnd = 1,2,3,4,7 or 8,
c                 but instead use mbdcnd = 5, 6, or 9.
c
c             2.  if b = pi, do not use mbdcnd = 1,2,3,4,5 or 6,
c                 but instead use mbdcnd = 7, 8, or 9.
c
c             when a = 0  and/or b = pi the only meaningful boundary
c             condition is du/dtheta = 0.  (see d. greenspan, 'numerical
c             analysis of elliptic boundary value problems,' harper and
c             row, 1965, chapter 5.)
c
c   bda
c     a one-dimensional array of length n that specifies the boundary
c     values (if any) of the solution at theta = a.  when
c     mbdcnd = 1, 2, or 7,
c
c              bda(j) = u(a,r(j)) ,              j=1,2,...,n.
c
c     when mbdcnd = 3, 4, or 8,
c
c              bda(j) = (d/dtheta)u(a,r(j)) ,    j=1,2,...,n.
c
c     when mbdcnd has any other value, bda is a dummy variable.
c
c   bdb
c     a one-dimensional array of length n that specifies the boundary
c     values of the solution at theta = b.  when mbdcnd = 1, 4, or 5,
c
c              bdb(j) = u(b,r(j)) ,              j=1,2,...,n.
c
c     when mbdcnd = 2,3, or 6,
c
c              bdb(j) = (d/dtheta)u(b,r(j)) ,    j=1,2,...,n.
c
c     when mbdcnd has any other value, bdb is a dummy variable.
c
c   c,d
c     the range of r , i.e. c .le. r .le. d.
c     c must be less than d.  c must be non-negative.
c
c   n
c     the number of unknowns in the interval (c,d).  the unknowns in
c     the r-direction are given by r(j) = c + (j-0.5)dr,
c     j=1,2,...,n, where dr = (d-c)/n.  n must be greater than 4.
c
c   nbdcnd
c     indicates the type of boundary conditions at r = c
c     and r = d.
c
c     = 1  if the solution is specified at r = c and r = d.
c
c     = 2  if the solution is specified at r = c and the derivative
c          of the solution with respect to r is specified at
c          r = d. (see note 1 below)
c
c     = 3  if the derivative of the solution with respect to r is
c          specified at r = c and r = d.
c
c     = 4  if the derivative of the solution with respect to r is
c          specified at r = c and the solution is specified at
c          r = d.
c
c     = 5  if the solution is unspecified at r = c = 0 (see note 2
c          below) and the solution is specified at r = d.
c
c     = 6  if the solution is unspecified at r = c = 0 (see note 2
c          below) and the derivative of the solution with respect to r
c          is specified at r = d.
c
c     note 1:  if c = 0 and mbdcnd = 3,6,8 or 9, the system of equations
c              to be solved is singular.  the unique solution is
c              determined by extrapolation to the specification of
c              u(theta(1),c).  but in these cases the right side of the
c              system will be perturbed by the constant pertrb.
c
c     note 2:  nbdcnd = 5 or 6 cannot be used with mbdcnd = 1, 2, 4, 5,
c              or 7 (the former indicates that the solution is
c              unspecified at r = 0; the latter indicates that the
c              solution is specified).  use instead nbdcnd = 1 or 2.
c
c   bdc
c     a one dimensional array of length m that specifies the boundary
c     values of the solution at r = c.   when nbdcnd = 1 or 2,
c
c              bdc(i) = u(theta(i),c) ,              i=1,2,...,m.
c
c     when nbdcnd = 3 or 4,
c
c              bdc(i) = (d/dr)u(theta(i),c),         i=1,2,...,m.
c
c     when nbdcnd has any other value, bdc is a dummy variable.
c
c   bdd
c     a one-dimensional array of length m that specifies the boundary
c     values of the solution at r = d.  when nbdcnd = 1 or 4,
c
c              bdd(i) = u(theta(i),d) ,              i=1,2,...,m.
c
c     when nbdcnd = 2 or 3,
c
c              bdd(i) = (d/dr)u(theta(i),d) ,        i=1,2,...,m.
c
c     when nbdcnd has any other value, bdd is a dummy variable.
c
c   elmbda
c     the constant lambda in the modified helmholtz equation.  if
c     lambda is greater than 0, a solution may not exist.  however,
c     hstcsp will attempt to find a solution.
c
c   f
c     a two-dimensional array that specifies the values of the right
c     side of the modified helmholtz equation.  for i=1,2,...,m and
c     j=1,2,...,n
c
c              f(i,j) = f(theta(i),r(j)) .
c
c     f must be dimensioned at least m x n.
c
c   idimf
c     the row (or first) dimension of the array f as it appears in the
c     program calling hstcsp.  this parameter is used to specify the
c     variable dimension of f.  idimf must be at least m.
c
c   w
c     a one-dimensional array that must be provided by the user for
c     work space.  with k = int(log2(n))+1 and l = 2**(k+1), w may
c     require up to (k-2)*l+k+max(2n,6m)+4(n+m)+5 locations.  the
c     actual number of locations used is computed by hstcsp and is
c     returned in the location w(1).
c
c
c            * * * * * *   on output   * * * * * *
c
c   f
c     contains the solution u(i,j) of the finite difference
c     approximation for the grid point (theta(i),r(j)) for
c     i=1,2,...,m, j=1,2,...,n.
c
c   pertrb
c     if a combination of periodic, derivative, or unspecified
c     boundary conditions is specified for a poisson equation
c     (lambda = 0), a solution may not exist.  pertrb is a con-
c     stant, calculated and subtracted from f, which ensures
c     that a solution exists.  hstcsp then computes this
c     solution, which is a least squares solution to the
c     original approximation.  this solution plus any constant is also
c     a solution; hence, the solution is not unique.  the value of
c     pertrb should be small compared to the right side f.
c     otherwise, a solution is obtained to an essentially different
c     problem.  this comparison should always be made to insure that
c     a meaningful solution has been obtained.
c
c   ierror
c     an error flag that indicates invalid input parameters.
c     except for numbers 0 and 10, a solution is not attempted.
c
c     =  0  no error
c
c     =  1  a .lt. 0 or b .gt. pi
c
c     =  2  a .ge. b
c
c     =  3  mbdcnd .lt. 1 or mbdcnd .gt. 9
c
c     =  4  c .lt. 0
c
c     =  5  c .ge. d
c
c     =  6  nbdcnd .lt. 1 or nbdcnd .gt. 6
c
c     =  7  n .lt. 5
c
c     =  8  nbdcnd = 5 or 6 and mbdcnd = 1, 2, 4, 5, or 7
c
c     =  9  c .gt. 0 and nbdcnd .ge. 5
c
c     = 10  elmbda .gt. 0
c
c     = 11  idimf .lt. m
c
c     = 12  m .lt. 5
c
c     = 13  a = 0 and mbdcnd =1,2,3,4,7 or 8
c
c     = 14  b = pi and mbdcnd .le. 6
c
c     = 15  a .gt. 0 and mbdcnd = 5, 6, or 9
c
c     = 16  b .lt. pi and mbdcnd .ge. 7
c
c     = 17  lambda .ne. 0 and nbdcnd .ge. 5
c
c     since this is the only means of indicating a possibly
c     incorrect call to hstcsp, the user should test ierror after
c     the call.
c
c   w
c     w(1) contains the required length of w.  also  w contains
c     intermediate values that must not be destroyed if hstcsp
c     will be called again with intl = 1.
c
c *long description:
c
c    * * * * * * *   program specifications    * * * * * * * * * * * *
c
c    dimension of   bda(n),bdb(n),bdc(m),bdd(m),f(idimf,n),
c    arguments      w(see argument list)
c
c    latest         june 1979
c    revision
c
c    subprograms    hstcsp,hstcs1,blktri,blktr1,indxa,indxb,indxc,
c    required       prod,prodp,cprod,cprodp,ppadd,psgf,bsrh,ppsgf,
c                   ppspf,compb,tevls,r1mach
c
c    special        none
c    conditions
c
c    common         cblkt
c    blocks
c
c    i/o            none
c
c    precision      single
c
c    specialist     roland sweet
c
c    language       fortran
c
c    history        written by roland sweet at ncar in may, 1977
c
c    algorithm      this subroutine defines the finite-difference
c                   equations, incorporates boundary data, adjusts the
c                   right side when the system is singular and calls
c                   blktri which solves the linear system of equations.
c
c    space          5269(decimal) = 12225(octal) locations on the
c    required       ncar control data 7600
c
c    timing and        the execution time t on the ncar control data
c    accuracy       7600 for subroutine hstcsp is roughly proportional
c                   to m*n*log2(n), but depends on the input parameter
c                   intl.  some values are listed in the table below.
c                      the solution process employed results in a loss
c                   of no more than four significant digits for n and m
c                   as large as 64.  more detailed information about
c                   accuracy can be found in the documentation for
c                   subroutine blktri which is the routine that
c                   actually solves the finite difference equations.
c
c
c                      m(=n)     intl      mbdcnd(=nbdcnd)     t(msecs)
c                      -----     ----      ---------------     --------
c
c                       32        0              1-6             132
c                       32        1              1-6              88
c                       64        0              1-6             546
c                       64        1              1-6             380
c
c    portability    american national standards institute fortran.
c                   the machine accuracy is set using function r1mach.
c
c    required       cos,sin,abs,sqrt
c    resident
c    routines
c
c    reference      swarztrauber, p.n., 'a direct method for the
c                   discrete solution of separable elliptic equations,'
c                   siam j. numer. anal. 11(1974), pp. 1136-1150.
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c***references  p. n. swarztrauber and r. sweet, efficient fortran
c                 subprograms for the solution of elliptic equations,
c                 ncar tn/ia-109, july 1975, 138 pp.
c               p. n. swarztrauber, a direct method for the discrete
c                 solution of separable elliptic equations, siam journal
c                 on numerical analysis 11, (1974), pp. 1136-1150.
c***routines called  hstcs1, pimach
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  hstcsp
c
c
      dimension       f(idimf,*) ,bda(*)     ,bdb(*)     ,bdc(*)     ,
     1                bdd(*)     ,w(*)
c***first executable statement  hstcsp
      pi = pimach(dum)
c
c     check for invalid input parameters
c
      ierror = 0
      if (a.lt.0. .or. b.gt.pi) ierror = 1
      if (a .ge. b) ierror = 2
      if (mbdcnd.lt.1 .or. mbdcnd.gt.9) ierror = 3
      if (c .lt. 0.) ierror = 4
      if (c .ge. d) ierror = 5
      if (nbdcnd.lt.1 .or. nbdcnd.gt.6) ierror = 6
      if (n .lt. 5) ierror = 7
      if ((nbdcnd.eq.5 .or. nbdcnd.eq.6) .and. (mbdcnd.eq.1 .or.
     1    mbdcnd.eq.2 .or. mbdcnd.eq.4 .or. mbdcnd.eq.5 .or.
     2                                                     mbdcnd.eq.7))
     3    ierror = 8
      if (c.gt.0. .and. nbdcnd.ge.5) ierror = 9
      if (idimf .lt. m) ierror = 11
      if (m .lt. 5) ierror = 12
      if (a.eq.0. .and. mbdcnd.ne.5 .and. mbdcnd.ne.6 .and. mbdcnd.ne.9)
     1    ierror = 13
      if (b.eq.pi .and. mbdcnd.le.6) ierror = 14
      if (a.gt.0. .and. (mbdcnd.eq.5 .or. mbdcnd.eq.6 .or. mbdcnd.eq.9))
     1    ierror = 15
      if (b.lt.pi .and. mbdcnd.ge.7) ierror = 16
      if (elmbda.ne.0. .and. nbdcnd.ge.5) ierror = 17
      if (ierror .ne. 0) go to 101
      iwbm = m+1
      iwcm = iwbm+m
      iwan = iwcm+m
      iwbn = iwan+n
      iwcn = iwbn+n
      iwsnth = iwcn+n
      iwrsq = iwsnth+m
      iwwrk = iwrsq+n
      ierr1 = 0
      call hstcs1 (intl,a,b,m,mbdcnd,bda,bdb,c,d,n,nbdcnd,bdc,bdd,
     1             elmbda,f,idimf,pertrb,ierr1,w,w(iwbm),w(iwcm),
     2             w(iwan),w(iwbn),w(iwcn),w(iwsnth),w(iwrsq),w(iwwrk))
      w(1) = w(iwwrk)+iwwrk-1
      ierror = ierr1
  101 continue
      return
      end
