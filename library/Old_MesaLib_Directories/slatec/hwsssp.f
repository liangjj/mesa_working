*deck hwsssp
      subroutine hwsssp (ts, tf, m, mbdcnd, bdts, bdtf, ps, pf, n,
     +   nbdcnd, bdps, bdpf, elmbda, f, idimf, pertrb, ierror, w)
c***begin prologue  hwsssp
c***purpose  solve a finite difference approximation to the helmholtz
c            equation in spherical coordinates and on the surface of the
c            unit sphere (radius of 1).
c***library   slatec (fishpack)
c***category  i2b1a1a
c***type      single precision (hwsssp-s)
c***keywords  elliptic, fishpack, helmholtz, pde, spherical
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c     subroutine hwsssp solves a finite difference approximation to the
c     helmholtz equation in spherical coordinates and on the surface of
c     the unit sphere (radius of 1):
c
c          (1/sin(theta))(d/dtheta)(sin(theta)(du/dtheta))
c
c             + (1/sin(theta)**2)(d/dphi)(du/dphi)
c
c             + lambda*u = f(theta,phi)
c
c     where theta is colatitude and phi is longitude.
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c     ts,tf
c       the range of theta (colatitude), i.e., ts .le. theta .le. tf.
c       ts must be less than tf.  ts and tf are in radians.  a ts of
c       zero corresponds to the north pole and a tf of pi corresponds to
c       the south pole.
c
c     * * * * * * * * * * * * * * important * * * * * * * * * * * * * *
c
c     if tf is equal to pi then it must be computed using the statement
c     tf = pimach(dum). this insures that tf in the users program is
c     equal to pi in this program which permits several tests of the
c     input parameters that otherwise would not be possible.
c
c
c     m
c       the number of panels into which the interval (ts,tf) is
c       subdivided.  hence, there will be m+1 grid points in the
c       theta-direction given by theta(i) = (i-1)dtheta+ts for
c       i = 1,2,...,m+1, where dtheta = (tf-ts)/m is the panel width.
c       m must be greater than 5.
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
c            specified at theta = ts (see note 1 below) and the
c            solution is unspecified at theta = tf = pi.
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
c            bdts(j) = (d/dtheta)u(ts,phi(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bdts is a dummy variable.
c
c     bdtf
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to theta at
c       theta = tf.  when mbdcnd = 2,3, or 6,
c
c            bdtf(j) = (d/dtheta)u(tf,phi(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bdtf is a dummy variable.
c
c     ps,pf
c       the range of phi (longitude), i.e., ps .le. phi .le. pf.  ps
c       must be less than pf.  ps and pf are in radians.  if ps = 0 and
c       pf = 2*pi, periodic boundary conditions are usually prescribed.
c
c     * * * * * * * * * * * * * * important * * * * * * * * * * * * * *
c
c     if pf is equal to 2*pi then it must be computed using the
c     statement pf = 2.*pimach(dum). this insures that pf in the users
c     program is equal to 2*pi in this program which permits tests of
c     the input parameters that otherwise would not be possible.
c
c
c     n
c       the number of panels into which the interval (ps,pf) is
c       subdivided.  hence, there will be n+1 grid points in the
c       phi-direction given by phi(j) = (j-1)dphi+ps  for
c       j = 1,2,...,n+1, where dphi = (pf-ps)/n is the panel width.
c       n must be greater than 4.
c
c     nbdcnd
c       indicates the type of boundary condition at phi = ps and
c       phi = pf.
c
c       = 0  if the solution is periodic in phi, i.e.,
c            u(i,j) = u(i,n+j).
c       = 1  if the solution is specified at phi = ps and phi = pf
c            (see note below).
c       = 2  if the solution is specified at phi = ps (see note below)
c            and the derivative of the solution with respect to phi is
c            specified at phi = pf.
c       = 3  if the derivative of the solution with respect to phi is
c            specified at phi = ps and phi = pf.
c       = 4  if the derivative of the solution with respect to phi is
c            specified at ps and the solution is specified at phi = pf
c            (see note below).
c
c       note:  nbdcnd = 1,2, or 4 cannot be used with
c              mbdcnd = 5,6,7,8, or 9 (the former indicates that the
c                       solution is specified at a pole, the latter
c                       indicates that the solution is unspecified).
c                       use instead
c              mbdcnd = 1 or 2  .
c
c     bdps
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to phi at
c       phi = ps.  when nbdcnd = 3 or 4,
c
c            bdps(i) = (d/dphi)u(theta(i),ps), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdps is a dummy variable.
c
c     bdpf
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to phi at
c       phi = pf.  when nbdcnd = 2 or 3,
c
c            bdpf(i) = (d/dphi)u(theta(i),pf), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdpf is a dummy variable.
c
c     elmbda
c       the constant lambda in the helmholtz equation.  if
c       lambda .gt. 0, a solution may not exist.  however, hwsssp will
c       attempt to find a solution.
c
c     f
c       a two-dimensional array that specifies the value of the right
c       side of the helmholtz equation and boundary values (if any).
c       for i = 2,3,...,m  and  j = 2,3,...,n
c
c            f(i,j) = f(theta(i),phi(j)).
c
c       on the boundaries f is defined by
c
c            mbdcnd   f(1,j)            f(m+1,j)
c            ------   ------------      ------------
c
c              1      u(ts,phi(j))      u(tf,phi(j))
c              2      u(ts,phi(j))      f(tf,phi(j))
c              3      f(ts,phi(j))      f(tf,phi(j))
c              4      f(ts,phi(j))      u(tf,phi(j))
c              5      f(0,ps)           u(tf,phi(j))   j = 1,2,...,n+1
c              6      f(0,ps)           f(tf,phi(j))
c              7      u(ts,phi(j))      f(pi,ps)
c              8      f(ts,phi(j))      f(pi,ps)
c              9      f(0,ps)           f(pi,ps)
c
c            nbdcnd   f(i,1)            f(i,n+1)
c            ------   --------------    --------------
c
c              0      f(theta(i),ps)    f(theta(i),ps)
c              1      u(theta(i),ps)    u(theta(i),pf)
c              2      u(theta(i),ps)    f(theta(i),pf)   i = 1,2,...,m+1
c              3      f(theta(i),ps)    f(theta(i),pf)
c              4      f(theta(i),ps)    u(theta(i),pf)
c
c       f must be dimensioned at least (m+1)*(n+1).
c
c      *note*
c
c       if the table calls for both the solution u and the right side f
c       at a corner then the solution must be specified.
c
c
c     idimf
c       the row (or first) dimension of the array f as it appears in the
c       program calling hwsssp.  this parameter is used to specify the
c       variable dimension of f.  idimf must be at least m+1  .
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space. w may require up to 4*(n+1)+(16+int(log2(n+1)))(m+1)
c       locations. the actual number of locations used is computed by
c       hwsssp and is output in location w(1). int( ) denotes the
c       fortran integer function.
c
c
c     * * * * * * * * * *     on output     * * * * * * * * * *
c
c     f
c       contains the solution u(i,j) of the finite difference
c       approximation for the grid point (theta(i),phi(j)),
c       i = 1,2,...,m+1,   j = 1,2,...,n+1  .
c
c     pertrb
c       if one specifies a combination of periodic, derivative or
c       unspecified boundary conditions for a poisson equation
c       (lambda = 0), a solution may not exist.  pertrb is a constant,
c       calculated and subtracted from f, which ensures that a solution
c       exists.  hwsssp then computes this solution, which is a least
c       squares solution to the original approximation.  this solution
c       is not unique and is unnormalized. the value of pertrb should
c       be small compared to the right side f. otherwise , a solution
c       is obtained to an essentially different problem. this comparison
c       should always be made to insure that a meaningful solution has
c       been obtained.
c
c     ierror
c       an error flag that indicates invalid input parameters.  except
c       for numbers 0 and 8, a solution is not attempted.
c
c       = 0  no error
c       = 1  ts.lt.0 or tf.gt.pi
c       = 2  ts.ge.tf
c       = 3  mbdcnd.lt.1 or mbdcnd.gt.9
c       = 4  ps.lt.0 or ps.gt.pi+pi
c       = 5  ps.ge.pf
c       = 6  n.lt.5
c       = 7  m.lt.5
c       = 8  nbdcnd.lt.0 or nbdcnd.gt.4
c       = 9  elmbda.gt.0
c       = 10 idimf.lt.m+1
c       = 11 nbdcnd equals 1,2 or 4 and mbdcnd.ge.5
c       = 12 ts.eq.0 and mbdcnd equals 3,4 or 8
c       = 13 tf.eq.pi and mbdcnd equals 2,3 or 6
c       = 14 mbdcnd equals 5,6 or 9 and ts.ne.0
c       = 15 mbdcnd.ge.7 and tf.ne.pi
c
c       since this is the only means of indicating a possibly incorrect
c       call to hwsssp, the user should test ierror after a call.
c
c     w
c       contains intermediate values that must not be destroyed if
c       hwsssp will be called again with intl = 1. w(1) contains the
c       required length of w .
c
c *long description:
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   bdts(n+1),bdtf(n+1),bdps(m+1),bdpf(m+1),
c     arguments      f(idimf,n+1),w(see argument list)
c
c     latest         january 1978
c     revision
c
c
c     subprograms    hwsssp,hwsss1,genbun,poisd2,poisn2,poisp2,cosgen,me
c     required       trix,tri3,pimach
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     paul swarztrauber
c
c     language       fortran
c
c     history        version 1 - september 1973
c                    version 2 - april     1976
c                    version 3 - january   1978
c
c     algorithm      the routine defines the finite difference
c                    equations, incorporates boundary data, and adjusts
c                    the right side of singular systems and then calls
c                    genbun to solve the system.
c
c     space
c     required       control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hwsssp is roughly proportional
c                    to m*n*log2(n), but also depends on the input
c                    parameters nbdcnd and mbdcnd.  some typical values
c                    are listed in the table below.
c                       the solution process employed results in a loss
c                    of no more than three significant digits for n and
c                    m as large as 64.  more detailed information about
c                    accuracy can be found in the documentation for
c                    subroutine genbun which is the routine that
c                    solves the finite difference equations.
c
c
c                       m(=n)    mbdcnd    nbdcnd    t(msecs)
c                       -----    ------    ------    --------
c
c                        32        0         0          31
c                        32        1         1          23
c                        32        3         3          36
c                        64        0         0         128
c                        64        1         1          96
c                        64        3         3         142
c
c     portability    american national standards institute fortran.
c                    the machine dependent constant pi is defined in
c                    function pimach.
c
c     required       sin,cos
c     resident
c     routines
c
c     references     p. n. swarztrauber,'the direct solution of the
c                    discrete poisson equation on the surface of a
c                    sphere, siam j. numer. anal.,15(1974), pp 212-215
c
c                    swarztrauber,p. and r. sweet, 'efficient fortran
c                    subprograms for the solution of elliptic equations'
c                    ncar tn/ia-109, july, 1975, 138 pp.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c***references  p. n. swarztrauber and r. sweet, efficient fortran
c                 subprograms for the solution of elliptic equations,
c                 ncar tn/ia-109, july 1975, 138 pp.
c               p. n. swarztrauber, the direct solution of the discrete
c                 poisson equation on the surface of a sphere, siam
c                 journal on numerical analysis 15 (1974), pp. 212-215.
c***routines called  hwsss1, pimach
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891009  removed unreferenced variable.  (wrb)
c   891009  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  hwsssp
c
      dimension       f(idimf,*) ,bdts(*)    ,bdtf(*)    ,bdps(*)    ,
     1                bdpf(*)    ,w(*)
c***first executable statement  hwsssp
      pi = pimach(dum)
      tpi = 2.*pi
      ierror = 0
      if (ts.lt.0. .or. tf.gt.pi) ierror = 1
      if (ts .ge. tf) ierror = 2
      if (mbdcnd.lt.1 .or. mbdcnd.gt.9) ierror = 3
      if (ps.lt.0. .or. pf.gt.tpi) ierror = 4
      if (ps .ge. pf) ierror = 5
      if (n .lt. 5) ierror = 6
      if (m .lt. 5) ierror = 7
      if (nbdcnd.lt.0 .or. nbdcnd.gt.4) ierror = 8
      if (elmbda .gt. 0.) ierror = 9
      if (idimf .lt. m+1) ierror = 10
      if ((nbdcnd.eq.1 .or. nbdcnd.eq.2 .or. nbdcnd.eq.4) .and.
     1    mbdcnd.ge.5) ierror = 11
      if (ts.eq.0. .and.
     1    (mbdcnd.eq.3 .or. mbdcnd.eq.4 .or. mbdcnd.eq.8)) ierror = 12
      if (tf.eq.pi .and.
     1    (mbdcnd.eq.2 .or. mbdcnd.eq.3 .or. mbdcnd.eq.6)) ierror = 13
      if ((mbdcnd.eq.5 .or. mbdcnd.eq.6 .or. mbdcnd.eq.9) .and.
     1    ts.ne.0.) ierror = 14
      if (mbdcnd.ge.7 .and. tf.ne.pi) ierror = 15
      if (ierror.ne.0 .and. ierror.ne.9) return
      call hwsss1 (ts,tf,m,mbdcnd,bdts,bdtf,ps,pf,n,nbdcnd,bdps,bdpf,
     1             elmbda,f,idimf,pertrb,w,w(m+2),w(2*m+3),w(3*m+4),
     2             w(4*m+5),w(5*m+6),w(6*m+7))
      w(1) = w(6*m+7)+6*(m+1)
      return
      end
