*deck hstssp
      subroutine hstssp (a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd,
     +   bdc, bdd, elmbda, f, idimf, pertrb, ierror, w)
c***begin prologue  hstssp
c***purpose  solve the standard five-point finite difference
c            approximation on a staggered grid to the helmholtz
c            equation in spherical coordinates and on the surface of
c            the unit sphere (radius of 1).
c***library   slatec (fishpack)
c***category  i2b1a1a
c***type      single precision (hstssp-s)
c***keywords  elliptic, fishpack, helmholtz, pde, spherical
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c     hstssp solves the standard five-point finite difference
c     approximation on a staggered grid to the helmholtz equation in
c     spherical coordinates and on the surface of the unit sphere
c     (radius of 1)
c
c             (1/sin(theta))(d/dtheta)(sin(theta)(du/dtheta)) +
c
c       (1/sin(theta)**2)(d/dphi)(du/dphi) + lambda*u = f(theta,phi)
c
c     where theta is colatitude and phi is longitude.
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    * * * * * * * *    parameter description     * * * * * * * * * *
c
c            * * * * * *   on input    * * * * * *
c
c   a,b
c     the range of theta (colatitude), i.e. a .le. theta .le. b.  a
c     must be less than b and a must be non-negative.  a and b are in
c     radians.  a = 0 corresponds to the north pole and b = pi
c     corresponds to the south pole.
c
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
c
c
c   m
c     the number of grid points in the interval (a,b).  the grid points
c     in the theta-direction are given by theta(i) = a + (i-0.5)dtheta
c     for i=1,2,...,m where dtheta =(b-a)/m.  m must be greater than 2.
c
c   mbdcnd
c     indicates the type of boundary conditions at theta = a and
c     theta = b.
c
c     = 1  if the solution is specified at theta = a and theta = b.
c          (see note 3 below)
c
c     = 2  if the solution is specified at theta = a and the derivative
c          of the solution with respect to theta is specified at
c          theta = b (see notes 2 and 3 below).
c
c     = 3  if the derivative of the solution with respect to theta is
c          specified at theta = a (see notes 1, 2 below) and theta = b.
c
c     = 4  if the derivative of the solution with respect to theta is
c          specified at theta = a (see notes 1 and 2 below) and the
c          solution is specified at theta = b.
c
c     = 5  if the solution is unspecified at theta = a = 0 and the
c          solution is specified at theta = b.  (see note 3 below)
c
c     = 6  if the solution is unspecified at theta = a = 0 and the
c          derivative of the solution with respect to theta is
c          specified at theta = b (see note 2 below).
c
c     = 7  if the solution is specified at theta = a and the
c          solution is unspecified at theta = b = pi. (see note 3 below)
c
c     = 8  if the derivative of the solution with respect to
c          theta is specified at theta = a (see note 1 below)
c          and the solution is unspecified at theta = b = pi.
c
c     = 9  if the solution is unspecified at theta = a = 0 and
c          theta = b = pi.
c
c     notes:  1.  if a = 0, do not use mbdcnd = 3, 4, or 8,
c                 but instead use mbdcnd = 5, 6, or 9.
c
c             2.  if b = pi, do not use mbdcnd = 2, 3, or 6,
c                 but instead use mbdcnd = 7, 8, or 9.
c
c             3.  when the solution is specified at theta = 0 and/or
c                 theta = pi and the other boundary conditions are
c                 combinations of unspecified, normal derivative, or
c                 periodicity a singular system results.  the unique
c                 solution is determined by extrapolation to the
c                 specification of the solution at either theta = 0 or
c                 theta = pi.  but in these cases the right side of the
c                 system will be perturbed by the constant pertrb.
c
c   bda
c     a one-dimensional array of length n that specifies the boundary
c     values (if any) of the solution at theta = a.  when
c     mbdcnd = 1, 2, or 7,
c
c              bda(j) = u(a,phi(j)) ,              j=1,2,...,n.
c
c     when mbdcnd = 3, 4, or 8,
c
c              bda(j) = (d/dtheta)u(a,phi(j)) ,    j=1,2,...,n.
c
c     when mbdcnd has any other value, bda is a dummy variable.
c
c   bdb
c     a one-dimensional array of length n that specifies the boundary
c     values of the solution at theta = b.  when mbdcnd = 1,4, or 5,
c
c              bdb(j) = u(b,phi(j)) ,              j=1,2,...,n.
c
c     when mbdcnd = 2,3, or 6,
c
c              bdb(j) = (d/dtheta)u(b,phi(j)) ,    j=1,2,...,n.
c
c     when mbdcnd has any other value, bdb is a dummy variable.
c
c   c,d
c     the range of phi (longitude), i.e. c .le. phi .le. d.
c     c must be less than d.  if d-c = 2*pi, periodic boundary
c     conditions are usually prescribed.
c
c   n
c     the number of unknowns in the interval (c,d).  the unknowns in
c     the phi-direction are given by phi(j) = c + (j-0.5)dphi,
c     j=1,2,...,n, where dphi = (d-c)/n.  n must be greater than 2.
c
c   nbdcnd
c     indicates the type of boundary conditions at phi = c
c     and phi = d.
c
c     = 0  if the solution is periodic in phi, i.e.
c          u(i,j) = u(i,n+j).
c
c     = 1  if the solution is specified at phi = c and phi = d
c          (see note below).
c
c     = 2  if the solution is specified at phi = c and the derivative
c          of the solution with respect to phi is specified at
c          phi = d (see note below).
c
c     = 3  if the derivative of the solution with respect to phi is
c          specified at phi = c and phi = d.
c
c     = 4  if the derivative of the solution with respect to phi is
c          specified at phi = c and the solution is specified at
c          phi = d (see note below).
c
c     note:  when nbdcnd = 1, 2, or 4, do not use mbdcnd = 5, 6, 7, 8,
c     or 9 (the former indicates that the solution is specified at
c     a pole; the latter indicates the solution is unspecified).  use
c     instead mbdcnd = 1 or 2.
c
c   bdc
c     a one dimensional array of length m that specifies the boundary
c     values of the solution at phi = c.   when nbdcnd = 1 or 2,
c
c              bdc(i) = u(theta(i),c) ,              i=1,2,...,m.
c
c     when nbdcnd = 3 or 4,
c
c              bdc(i) = (d/dphi)u(theta(i),c),       i=1,2,...,m.
c
c     when nbdcnd = 0, bdc is a dummy variable.
c
c   bdd
c     a one-dimensional array of length m that specifies the boundary
c     values of the solution at phi = d.  when nbdcnd = 1 or 4,
c
c              bdd(i) = u(theta(i),d) ,              i=1,2,...,m.
c
c     when nbdcnd = 2 or 3,
c
c              bdd(i) = (d/dphi)u(theta(i),d) ,      i=1,2,...,m.
c
c     when nbdcnd = 0, bdd is a dummy variable.
c
c   elmbda
c     the constant lambda in the helmholtz equation.  if lambda is
c     greater than 0, a solution may not exist.  however, hstssp will
c     attempt to find a solution.
c
c   f
c     a two-dimensional array that specifies the values of the right
c     side of the helmholtz equation.  for i=1,2,...,m and j=1,2,...,n
c
c              f(i,j) = f(theta(i),phi(j)) .
c
c     f must be dimensioned at least m x n.
c
c   idimf
c     the row (or first) dimension of the array f as it appears in the
c     program calling hstssp.  this parameter is used to specify the
c     variable dimension of f.  idimf must be at least m.
c
c   w
c     a one-dimensional array that must be provided by the user for
c     work space.  w may require up to 13m + 4n + m*int(log2(n))
c     locations.  the actual number of locations used is computed by
c     hstssp and is returned in the location w(1).
c
c
c            * * * * * *   on output   * * * * * *
c
c   f
c     contains the solution u(i,j) of the finite difference
c     approximation for the grid point (theta(i),phi(j)) for
c     i=1,2,...,m, j=1,2,...,n.
c
c   pertrb
c     if a combination of periodic, derivative, or unspecified
c     boundary conditions is specified for a poisson equation
c     (lambda = 0), a solution may not exist.  pertrb is a con-
c     stant, calculated and subtracted from f, which ensures
c     that a solution exists.  hstssp then computes this
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
c      except for numbers 0 and 14, a solution is not attempted.
c
c     =  0  no error
c
c     =  1  a .lt. 0 or b .gt. pi
c
c     =  2  a .ge. b
c
c     =  3  mbdcnd .lt. 1 or mbdcnd .gt. 9
c
c     =  4  c .ge. d
c
c     =  5  n .le. 2
c
c     =  6  nbdcnd .lt. 0 or nbdcnd .gt. 4
c
c     =  7  a .gt. 0 and mbdcnd = 5, 6, or 9
c
c     =  8  a = 0 and mbdcnd = 3, 4, or 8
c
c     =  9  b .lt. pi and mbdcnd .ge. 7
c
c     = 10  b = pi and mbdcnd = 2,3, or 6
c
c     = 11  mbdcnd .ge. 5 and ndbcnd = 1, 2, or 4
c
c     = 12  idimf .lt. m
c
c     = 13  m .le. 2
c
c     = 14  lambda .gt. 0
c
c     since this is the only means of indicating a possibly
c     incorrect call to hstssp, the user should test ierror after
c     the call.
c
c   w
c     w(1) contains the required length of w.
c
c *long description:
c
c    * * * * * * *   program specifications    * * * * * * * * * * * *
c
c    dimension of   bda(n),bdb(n),bdc(m),bdd(m),f(idimf,n),
c    arguments      w(see argument list)
c
c    latest         june 1, 1977
c    revision
c
c    subprograms    hstssp,poistg,postg2,genbun,poisd2,poisn2,poisp2,
c    required       cosgen,merge,trix,tri3,pimach
c
c    special        none
c    conditions
c
c    common         none
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
c    history        written by roland sweet at ncar in april, 1977
c
c    algorithm      this subroutine defines the finite-difference
c                   equations, incorporates boundary data, adjusts the
c                   right side when the system is singular and calls
c                   either poistg or genbun which solves the linear
c                   system of equations.
c
c    space          8427(decimal) = 20353(octal) locations on the
c    required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hstssp is roughly proportional
c                    to m*n*log2(n).  some typical values are listed in
c                    the table below.
c                       the solution process employed results in a loss
c                    of no more than four significant digits for n and m
c                    as large as 64.  more detailed information about
c                    accuracy can be found in the documentation for
c                    subroutine poistg which is the routine that
c                    actually solves the finite difference equations.
c
c
c                       m(=n)    mbdcnd    nbdcnd    t(msecs)
c                       -----    ------    ------    --------
c
c                        32       1-9       1-4         56
c                        64       1-9       1-4        230
c
c    portability     american national standards institute fortran.
c                    the machine dependent constant pi is defined in
c                    function pimach.
c
c    required       cos
c    resident
c    routines
c
c    reference      schumann, u. and r. sweet,'a direct method for
c                   the solution of poisson's equation with neumann
c                   boundary conditions on a staggered grid of
c                   arbitrary size,' j. comp. phys. 20(1976),
c                   pp. 171-182.
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c***references  u. schumann and r. sweet, a direct method for the
c                 solution of poisson's equation with neumann boundary
c                 conditions on a staggered grid of arbitrary size,
c                 journal of computational physics 20, (1976),
c                 pp. 171-182.
c***routines called  genbun, pimach, poistg
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  hstssp
c
c
      dimension       f(idimf,*) ,bda(*)     ,bdb(*)     ,bdc(*)     ,
     1                bdd(*)     ,w(*)
c***first executable statement  hstssp
      ierror = 0
      pi = pimach(dum)
      if (a.lt.0. .or. b.gt.pi) ierror = 1
      if (a .ge. b) ierror = 2
      if (mbdcnd.le.0 .or. mbdcnd.gt.9) ierror = 3
      if (c .ge. d) ierror = 4
      if (n .le. 2) ierror = 5
      if (nbdcnd.lt.0 .or. nbdcnd.ge.5) ierror = 6
      if (a.gt.0. .and. (mbdcnd.eq.5 .or. mbdcnd.eq.6 .or. mbdcnd.eq.9))
     1    ierror = 7
      if (a.eq.0. .and. (mbdcnd.eq.3 .or. mbdcnd.eq.4 .or. mbdcnd.eq.8))
     1    ierror = 8
      if (b.lt.pi .and. mbdcnd.ge.7) ierror = 9
      if (b.eq.pi .and. (mbdcnd.eq.2 .or. mbdcnd.eq.3 .or. mbdcnd.eq.6))
     1    ierror = 10
      if (mbdcnd.ge.5 .and.
     1    (nbdcnd.eq.1 .or. nbdcnd.eq.2 .or. nbdcnd.eq.4)) ierror = 11
      if (idimf .lt. m) ierror = 12
      if (m .le. 2) ierror = 13
      if (ierror .ne. 0) return
      deltar = (b-a)/m
      dlrsq = deltar**2
      deltht = (d-c)/n
      dlthsq = deltht**2
      np = nbdcnd+1
      isw = 1
      jsw = 1
      mb = mbdcnd
      if (elmbda .ne. 0.) go to 105
      go to (101,102,105,103,101,105,101,105,105),mbdcnd
  101 if (a.ne.0. .or. b.ne.pi) go to 105
      mb = 9
      go to 104
  102 if (a .ne. 0.) go to 105
      mb = 6
      go to 104
  103 if (b .ne. pi) go to 105
      mb = 8
  104 jsw = 2
  105 continue
c
c     define a,b,c coefficients in w-array.
c
      iwb = m
      iwc = iwb+m
      iwr = iwc+m
      iws = iwr+m
      do 106 i=1,m
         j = iwr+i
         w(j) = sin(a+(i-0.5)*deltar)
         w(i) = sin((a+(i-1)*deltar))/dlrsq
  106 continue
      mm1 = m-1
      do 107 i=1,mm1
         k = iwc+i
         w(k) = w(i+1)
         j = iwr+i
         k = iwb+i
         w(k) = elmbda*w(j)-(w(i)+w(i+1))
  107 continue
      w(iwr) = sin(b)/dlrsq
      w(iwc) = elmbda*w(iws)-(w(m)+w(iwr))
      do 109 i=1,m
         j = iwr+i
         a1 = w(j)
         do 108 j=1,n
            f(i,j) = a1*f(i,j)
  108    continue
  109 continue
c
c     enter boundary data for theta-boundaries.
c
      go to (110,110,112,112,114,114,110,112,114),mb
  110 a1 = 2.*w(1)
      w(iwb+1) = w(iwb+1)-w(1)
      do 111 j=1,n
         f(1,j) = f(1,j)-a1*bda(j)
  111 continue
      go to 114
  112 a1 = deltar*w(1)
      w(iwb+1) = w(iwb+1)+w(1)
      do 113 j=1,n
         f(1,j) = f(1,j)+a1*bda(j)
  113 continue
  114 go to (115,117,117,115,115,117,119,119,119),mb
  115 a1 = 2.*w(iwr)
      w(iwc) = w(iwc)-w(iwr)
      do 116 j=1,n
         f(m,j) = f(m,j)-a1*bdb(j)
  116 continue
      go to 119
  117 a1 = deltar*w(iwr)
      w(iwc) = w(iwc)+w(iwr)
      do 118 j=1,n
         f(m,j) = f(m,j)-a1*bdb(j)
  118 continue
c
c     enter boundary data for phi-boundaries.
c
  119 a1 = 2./dlthsq
      go to (129,120,120,122,122),np
  120 do 121 i=1,m
         j = iwr+i
         f(i,1) = f(i,1)-a1*bdc(i)/w(j)
  121 continue
      go to 124
  122 a1 = 1./deltht
      do 123 i=1,m
         j = iwr+i
         f(i,1) = f(i,1)+a1*bdc(i)/w(j)
  123 continue
  124 a1 = 2./dlthsq
      go to (129,125,127,127,125),np
  125 do 126 i=1,m
         j = iwr+i
         f(i,n) = f(i,n)-a1*bdd(i)/w(j)
  126 continue
      go to 129
  127 a1 = 1./deltht
      do 128 i=1,m
         j = iwr+i
         f(i,n) = f(i,n)-a1*bdd(i)/w(j)
  128 continue
  129 continue
c
c     adjust right side of singular problems to insure existence of a
c     solution.
c
      pertrb = 0.
      if (elmbda) 139,131,130
  130 ierror = 14
      go to 139
  131 go to (139,139,132,139,139,132,139,132,132),mb
  132 go to (133,139,139,133,139),np
  133 continue
      isw = 2
      do 135 j=1,n
         do 134 i=1,m
            pertrb = pertrb+f(i,j)
  134    continue
  135 continue
      a1 = n*(cos(a)-cos(b))/(2.*sin(0.5*deltar))
      pertrb = pertrb/a1
      do 137 i=1,m
         j = iwr+i
         a1 = pertrb*w(j)
         do 136 j=1,n
            f(i,j) = f(i,j)-a1
  136    continue
  137 continue
      a2 = 0.
      a3 = 0.
      do 138 j=1,n
         a2 = a2+f(1,j)
         a3 = a3+f(m,j)
  138 continue
      a2 = a2/w(iwr+1)
      a3 = a3/w(iws)
  139 continue
c
c     multiply i-th equation through by  r(i)*deltht**2
c
      do 141 i=1,m
         j = iwr+i
         a1 = dlthsq*w(j)
         w(i) = a1*w(i)
         j = iwc+i
         w(j) = a1*w(j)
         j = iwb+i
         w(j) = a1*w(j)
         do 140 j=1,n
            f(i,j) = a1*f(i,j)
  140    continue
  141 continue
      lp = nbdcnd
      w(1) = 0.
      w(iwr) = 0.
c
c     call poistg or genbun to solve the system of equations.
c
      if (nbdcnd .eq. 0) go to 142
      call poistg (lp,n,1,m,w,w(iwb+1),w(iwc+1),idimf,f,ierr1,w(iwr+1))
      go to 143
  142 call genbun (lp,n,1,m,w,w(iwb+1),w(iwc+1),idimf,f,ierr1,w(iwr+1))
  143 continue
      w(1) = w(iwr+1)+3*m
      if (isw.ne.2 .or. jsw.ne.2) go to 150
      if (mb .ne. 8) go to 145
      a1 = 0.
      do 144 j=1,n
         a1 = a1+f(m,j)
  144 continue
      a1 = (a1-dlrsq*a3/16.)/n
      if (nbdcnd .eq. 3) a1 = a1+(bdd(m)-bdc(m))/(d-c)
      a1 = bdb(1)-a1
      go to 147
  145 a1 = 0.
      do 146 j=1,n
         a1 = a1+f(1,j)
  146 continue
      a1 = (a1-dlrsq*a2/16.)/n
      if (nbdcnd .eq. 3) a1 = a1+(bdd(1)-bdc(1))/(d-c)
      a1 = bda(1)-a1
  147 do 149 i=1,m
         do 148 j=1,n
            f(i,j) = f(i,j)+a1
  148    continue
  149 continue
  150 continue
      return
      end
