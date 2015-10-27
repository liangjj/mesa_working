*deck hwscsp
      subroutine hwscsp (intl, ts, tf, m, mbdcnd, bdts, bdtf, rs, rf, n,
     +   nbdcnd, bdrs, bdrf, elmbda, f, idimf, pertrb, ierror, w)
c***begin prologue  hwscsp
c***purpose  solve a finite difference approximation to the modified
c            helmholtz equation in spherical coordinates assuming
c            axisymmetry  (no dependence on longitude).
c***library   slatec (fishpack)
c***category  i2b1a1a
c***type      single precision (hwscsp-s)
c***keywords  elliptic, fishpack, helmholtz, pde, spherical
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c     subroutine hwscsp solves a finite difference approximation to the
c       modified helmholtz equation in spherical coordinates assuming
c       axisymmetry  (no dependence on longitude)
c
c          (1/r**2)(d/dr)((r**2)(d/dr)u)
c
c             + (1/(r**2)sin(theta))(d/dtheta)(sin(theta)(d/dtheta)u)
c
c             + (lambda/(rsin(theta))**2)u = f(theta,r).
c
c     this two dimensional modified helmholtz equation results from
c     the fourier transform of the three dimensional poisson equation
c
c     * * * * * * * * * *     on input     * * * * * * * * * *
c
c     intl
c       = 0  on initial entry to hwscsp or if any of the arguments
c            rs, rf, n, nbdcnd are changed from a previous call.
c       = 1  if rs, rf, n, nbdcnd are all unchanged from previous call
c            to hwscsp.
c
c       note   a call with intl=0 takes approximately 1.5 times as
c              much time as a call with intl = 1.  once a call with
c              intl = 0 has been made then subsequent solutions
c              corresponding to different f, bdts, bdtf, bdrs, bdrf can
c              be obtained faster with intl = 1 since initialization is
c              not repeated.
c
c     ts,tf
c       the range of theta (colatitude), i.e., ts .le. theta .le. tf.
c       ts must be less than tf.  ts and tf are in radians.  a ts of
c       zero corresponds to the north pole and a tf of pi corresponds
c       to the south pole.
c
c     * * * * * * * * * * * * * * important * * * * * * * * * * * * * *
c
c     if tf is equal to pi then it must be computed using the statement
c     tf = pimach(dum). this insures that tf in the users program is
c     equal to pi in this program which permits several tests of the
c     input parameters that otherwise would not be possible.
c
c     m
c       the number of panels into which the interval (ts,tf) is
c       subdivided.  hence, there will be m+1 grid points in the
c       theta-direction given by theta(k) = (i-1)dtheta+ts for
c       i = 1,2,...,m+1, where dtheta = (tf-ts)/m is the panel width.
c
c     mbdcnd
c       indicates the type of boundary condition at theta = ts and
c       theta = tf.
c
c       = 1  if the solution is specified at theta = ts and theta = tf.
c       = 2  if the solution is specified at theta = ts and the
c            derivative of the solution with respect to theta is
c            specified at theta = tf (see note 2 below).
c       = 3  if the derivative of the solution with respect to theta is
c            specified at theta = ts and theta = tf (see notes 1,2
c            below).
c       = 4  if the derivative of the solution with respect to theta is
c            specified at theta = ts (see note 1 below) and the
c            solution is specified at theta = tf.
c       = 5  if the solution is unspecified at theta = ts = 0 and the
c            solution is specified at theta = tf.
c       = 6  if the solution is unspecified at theta = ts = 0 and the
c            derivative of the solution with respect to theta is
c            specified at theta = tf (see note 2 below).
c       = 7  if the solution is specified at theta = ts and the
c            solution is unspecified at theta = tf = pi.
c       = 8  if the derivative of the solution with respect to theta is
c            specified at theta = ts (see note 1 below) and the solution
c            is unspecified at theta = tf = pi.
c       = 9  if the solution is unspecified at theta = ts = 0 and
c            theta = tf = pi.
c
c       notes:  1.  if ts = 0, do not use mbdcnd = 3,4, or 8, but
c                   instead use mbdcnd = 5,6, or 9  .
c               2.  if tf = pi, do not use mbdcnd = 2,3, or 6, but
c                   instead use mbdcnd = 7,8, or 9  .
c
c     bdts
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to theta at
c       theta = ts.  when mbdcnd = 3,4, or 8,
c
c            bdts(j) = (d/dtheta)u(ts,r(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bdts is a dummy variable.
c
c     bdtf
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to theta at
c       theta = tf.  when mbdcnd = 2,3, or 6,
c
c            bdtf(j) = (d/dtheta)u(tf,r(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bdtf is a dummy variable.
c
c     rs,rf
c       the range of r, i.e., rs .le. r .lt. rf.  rs must be less than
c       rf.  rs must be non-negative.
c
c       n
c       the number of panels into which the interval (rs,rf) is
c       subdivided.  hence, there will be n+1 grid points in the
c       r-direction given by r(j) = (j-1)dr+rs for j = 1,2,...,n+1,
c       where dr = (rf-rs)/n is the panel width.
c       n must be greater than 2
c
c     nbdcnd
c       indicates the type of boundary condition at r = rs and r = rf.
c
c       = 1  if the solution is specified at r = rs and r = rf.
c       = 2  if the solution is specified at r = rs and the derivative
c            of the solution with respect to r is specified at r = rf.
c       = 3  if the derivative of the solution with respect to r is
c            specified at r = rs and r = rf.
c       = 4  if the derivative of the solution with respect to r is
c            specified at rs and the solution is specified at r = rf.
c       = 5  if the solution is unspecified at r = rs = 0 (see note
c            below) and the solution is specified at r = rf.
c       = 6  if the solution is unspecified at r = rs = 0 (see note
c            below) and the derivative of the solution with respect to
c            r is specified at r = rf.
c
c       note:  nbdcnd = 5 or 6 cannot be used with
c              mbdcnd = 1,2,4,5, or 7 (the former indicates that the
c                       solution is unspecified at r = 0, the latter
c                       indicates that the solution is specified).
c                       use instead
c              nbdcnd = 1 or 2  .
c
c     bdrs
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to r at r = rs.
c       when nbdcnd = 3 or 4,
c
c            bdrs(i) = (d/dr)u(theta(i),rs), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdrs is a dummy variable.
c
c     bdrf
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to r at r = rf.
c       when nbdcnd = 2,3, or 6,
c
c            bdrf(i) = (d/dr)u(theta(i),rf), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdrf is a dummy variable.
c
c     elmbda
c       the constant lambda in the helmholtz equation.  if
c       lambda .gt. 0, a solution may not exist.  however, hwscsp will
c       attempt to find a solution.  if nbdcnd = 5 or 6 or
c       mbdcnd = 5,6,7,8, or 9, elmbda must be zero.
c
c     f
c       a two-dimensional array that specifies the value of the right
c       side of the helmholtz equation and boundary values (if any).
c       for i = 2,3,...,m and j = 2,3,...,n
c
c            f(i,j) = f(theta(i),r(j)).
c
c       on the boundaries f is defined by
c
c            mbdcnd   f(1,j)            f(m+1,j)
c            ------   ----------        ----------
c
c              1      u(ts,r(j))        u(tf,r(j))
c              2      u(ts,r(j))        f(tf,r(j))
c              3      f(ts,r(j))        f(tf,r(j))
c              4      f(ts,r(j))        u(tf,r(j))
c              5      f(0,r(j))         u(tf,r(j))   j = 1,2,...,n+1
c              6      f(0,r(j))         f(tf,r(j))
c              7      u(ts,r(j))        f(pi,r(j))
c              8      f(ts,r(j))        f(pi,r(j))
c              9      f(0,r(j))         f(pi,r(j))
c
c            nbdcnd   f(i,1)            f(i,n+1)
c            ------   --------------    --------------
c
c              1      u(theta(i),rs)    u(theta(i),rf)
c              2      u(theta(i),rs)    f(theta(i),rf)
c              3      f(theta(i),rs)    f(theta(i),rf)
c              4      f(theta(i),rs)    u(theta(i),rf)   i = 1,2,...,m+1
c              5      f(ts,0)           u(theta(i),rf)
c              6      f(ts,0)           f(theta(i),rf)
c
c       f must be dimensioned at least (m+1)*(n+1).
c
c       note
c
c       if the table calls for both the solution u and the right side f
c       at a corner then the solution must be specified.
c
c     idimf
c       the row (or first) dimension of the array f as it appears in the
c       program calling hwscsp.  this parameter is used to specify the
c       variable dimension of f.  idimf must be at least m+1  .
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space. its length can be computed from the formula below
c       which depends on the value of nbdcnd.
c
c       if nbdcnd=2,4 or 6 define nunk=n
c       if nbdcnd=1 or 5   define nunk=n-1
c       if nbdcnd=3        define nunk=n+1
c
c       now set k=int(log2(nunk))+1 and l=2**(k+1) then w must be
c       dimensioned at least (k-2)*l+k+5*(m+n)+max(2*n,6*m)+23
c
c       **important** for purposes of checking, the required length
c                     of w is computed by hwscsp and stored in w(1)
c                     in floating point format.
c
c
c     * * * * * * * * * *     on output     * * * * * * * * * *
c
c     f
c       contains the solution u(i,j) of the finite difference
c       approximation for the grid point (theta(i),r(j)),
c       i = 1,2,...,m+1,   j = 1,2,...,n+1  .
c
c     pertrb
c       if a combination of periodic or derivative boundary conditions
c       is specified for a poisson equation (lambda = 0), a solution may
c       not exist.  pertrb is a constant, calculated and subtracted from
c       f, which ensures that a solution exists.  hwscsp then computes
c       this solution, which is a least squares solution to the original
c       approximation. this solution is not unique and is unnormalized.
c       the value of pertrb should be small compared to the right side
c       f. otherwise , a solution is obtained to an essentially
c       different problem. this comparison should always be made to
c       insure that a meaningful solution has been obtained.
c
c     ierror
c       an error flag that indicates invalid input parameters.  except
c       for numbers 0 and 10, a solution is not attempted.
c
c       = 1  ts.lt.0. or tf.gt.pi
c       = 2  ts.ge.tf
c       = 3  m.lt.5
c       = 4  mbdcnd.lt.1 or mbdcnd.gt.9
c       = 5  rs.lt.0
c       = 6  rs.ge.rf
c       = 7  n.lt.5
c       = 8  nbdcnd.lt.1 or nbdcnd.gt.6
c       = 9  elmbda.gt.0
c       = 10 idimf.lt.m+1
c       = 11 elmbda.ne.0 and mbdcnd.ge.5
c       = 12 elmbda.ne.0 and nbdcnd equals 5 or 6
c       = 13 mbdcnd equals 5,6 or 9 and ts.ne.0
c       = 14 mbdcnd.ge.7 and tf.ne.pi
c       = 15 ts.eq.0 and mbdcnd equals 3,4 or 8
c       = 16 tf.eq.pi and mbdcnd equals 2,3 or 6
c       = 17 nbdcnd.ge.5 and rs.ne.0
c       = 18 nbdcnd.ge.5 and mbdcnd equals 1,2,4,5 or 7
c
c       since this is the only means of indicating a possibly incorrect
c       call to hwscsp, the user should test ierror after a call.
c
c     w
c       contains intermediate values that must not be destroyed if
c       hwscsp will be called again with intl = 1.  w(1) contains the
c       number of locations which w must have.
c
c *long description:
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   bdts(n+1),bdtf(n+1),bdrs(m+1),bdrf(m+1),
c     arguments      f(idimf,n+1),w(see argument list)
c
c     latest         june 1979
c     revision
c
c     subprograms    hwscsp,hwscs1,blktri,blktr1,prod,prodp,cprod,cprodp
c     required       ,combp,ppadd,psgf,bsrh,ppsgf,ppspf,tevls,indxa,
c                    ,indxb,indxc,r1mach
c
c     special
c     conditions
c
c     common         cblkt
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     paul n swarztrauber
c
c     language       fortran
c
c     history        version 1 september 1973
c                    version 2 april     1976
c                    version 3 june      1979
c
c     algorithm      the routine defines the finite difference
c                    equations, incorporates boundary data, and adjusts
c                    the right side of singular systems and then calls
c                    blktri to solve the system.
c
c     space
c     required
c
c     portability    american national standards institute fortran.
c                    the machine accuracy is set using function r1mach.
c
c     required       none
c     resident
c     routines
c
c     reference      swarztrauber,p. and r. sweet, 'efficient fortran
c                    subprograms for the solution of elliptic equations'
c                    ncar tn/ia-109, july, 1975, 138 pp.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c***references  p. n. swarztrauber and r. sweet, efficient fortran
c                 subprograms for the solution of elliptic equations,
c                 ncar tn/ia-109, july 1975, 138 pp.
c***routines called  hwscs1, pimach
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  hwscsp
c
      dimension       f(idimf,*) ,bdts(*)    ,bdtf(*)    ,bdrs(*)    ,
     1                bdrf(*)    ,w(*)
c***first executable statement  hwscsp
      pi = pimach(dum)
      ierror = 0
      if (ts.lt.0. .or. tf.gt.pi) ierror = 1
      if (ts .ge. tf) ierror = 2
      if (m .lt. 5) ierror = 3
      if (mbdcnd.lt.1 .or. mbdcnd.gt.9) ierror = 4
      if (rs .lt. 0.) ierror = 5
      if (rs .ge. rf) ierror = 6
      if (n .lt. 5) ierror = 7
      if (nbdcnd.lt.1 .or. nbdcnd.gt.6) ierror = 8
      if (elmbda .gt. 0.) ierror = 9
      if (idimf .lt. m+1) ierror = 10
      if (elmbda.ne.0. .and. mbdcnd.ge.5) ierror = 11
      if (elmbda.ne.0. .and. (nbdcnd.eq.5 .or. nbdcnd.eq.6)) ierror = 12
      if ((mbdcnd.eq.5 .or. mbdcnd.eq.6 .or. mbdcnd.eq.9) .and.
     1    ts.ne.0.) ierror = 13
      if (mbdcnd.ge.7 .and. tf.ne.pi) ierror = 14
      if (ts.eq.0. .and.
     1    (mbdcnd.eq.4 .or. mbdcnd.eq.8 .or. mbdcnd.eq.3)) ierror = 15
      if (tf.eq.pi .and.
     1    (mbdcnd.eq.2 .or. mbdcnd.eq.3 .or. mbdcnd.eq.6)) ierror = 16
      if (nbdcnd.ge.5 .and. rs.ne.0.) ierror = 17
      if (nbdcnd.ge.5 .and. (mbdcnd.eq.1 .or. mbdcnd.eq.2 .or.
     1                                    mbdcnd.eq.5 .or. mbdcnd.eq.7))
     2    ierror = 18
      if (ierror.ne.0 .and. ierror.ne.9) return
      nck = n
      go to (101,103,102,103,101,103),nbdcnd
  101 nck = nck-1
      go to 103
  102 nck = nck+1
  103 l = 2
      k = 1
  104 l = l+l
      k = k+1
      if (nck-l) 105,105,104
  105 l = l+l
      np1 = n+1
      mp1 = m+1
      i1 = (k-2)*l+k+max(2*n,6*m)+13
      i2 = i1+np1
      i3 = i2+np1
      i4 = i3+np1
      i5 = i4+np1
      i6 = i5+np1
      i7 = i6+mp1
      i8 = i7+mp1
      i9 = i8+mp1
      i10 = i9+mp1
      w(1) = i10+m
      call hwscs1 (intl,ts,tf,m,mbdcnd,bdts,bdtf,rs,rf,n,nbdcnd,bdrs,
     1             bdrf,elmbda,f,idimf,pertrb,w(2),w(i1),w(i2),w(i3),
     2             w(i4),w(i5),w(i6),w(i7),w(i8),w(i9),w(i10))
      return
      end
