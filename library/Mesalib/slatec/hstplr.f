*deck hstplr
      subroutine hstplr (a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd,
     +   bdc, bdd, elmbda, f, idimf, pertrb, ierror, w)
c***begin prologue  hstplr
c***purpose  solve the standard five-point finite difference
c            approximation on a staggered grid to the helmholtz equation
c            in polar coordinates.
c***library   slatec (fishpack)
c***category  i2b1a1a
c***type      single precision (hstplr-s)
c***keywords  elliptic, fishpack, helmholtz, pde, polar
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c      hstplr solves the standard five-point finite difference
c      approximation on a staggered grid to the helmholtz equation in
c      polar coordinates
c
c      (1/r)(d/dr)(r(du/dr)) + (1/r**2)(d/dtheta)(du/dtheta)
c
c                      + lambda*u = f(r,theta)
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c    a,b
c      the range of r, i.e. a .le. r .le. b.  a must be less than b and
c      a must be non-negative.
c
c    m
c      the number of grid points in the interval (a,b).  the grid points
c      in the r-direction are given by r(i) = a + (i-0.5)dr for
c      i=1,2,...,m where dr =(b-a)/m.  m must be greater than 2.
c
c    mbdcnd
c      indicates the type of boundary conditions at r = a and r = b.
c
c      = 1  if the solution is specified at r = a and r = b.
c
c      = 2  if the solution is specified at r = a and the derivative
c           of the solution with respect to r is specified at r = b.
c           (see note 1 below)
c
c      = 3  if the derivative of the solution with respect to r is
c           specified at r = a (see note 2 below) and r = b.
c
c      = 4  if the derivative of the solution with respect to r is
c           specified at r = a (see note 2 below) and the solution is
c           specified at r = b.
c
c      = 5  if the solution is unspecified at r = a = 0 and the solution
c           is specified at r = b.
c
c      = 6  if the solution is unspecified at r = a = 0 and the
c           derivative of the solution with respect to r is specified at
c           r = b.
c
c      note 1:  if a = 0, mbdcnd = 2, and nbdcnd = 0 or 3, the system of
c               equations to be solved is singular.  the unique solution
c               is determined by extrapolation to the specification of
c               u(0,theta(1)).  but in this case the right side of the
c               system will be perturbed by the constant pertrb.
c
c      note 2:  if a = 0, do not use mbdcnd = 3 or 4, but instead use
c               mbdcnd = 1,2,5, or 6.
c
c    bda
c      a one-dimensional array of length n that specifies the boundary
c      values (if any) of the solution at r = a.  when mbdcnd = 1 or 2,
c
c               bda(j) = u(a,theta(j)) ,          j=1,2,...,n.
c
c      when mbdcnd = 3 or 4,
c
c               bda(j) = (d/dr)u(a,theta(j)) ,    j=1,2,...,n.
c
c      when mbdcnd = 5 or 6, bda is a dummy variable.
c
c    bdb
c      a one-dimensional array of length n that specifies the boundary
c      values of the solution at r = b.  when mbdcnd = 1,4, or 5,
c
c               bdb(j) = u(b,theta(j)) ,          j=1,2,...,n.
c
c      when mbdcnd = 2,3, or 6,
c
c               bdb(j) = (d/dr)u(b,theta(j)) ,    j=1,2,...,n.
c
c    c,d
c      the range of theta, i.e. c .le. theta .le. d.  c must be less
c      than d.
c
c    n
c      the number of unknowns in the interval (c,d).  the unknowns in
c      the theta-direction are given by theta(j) = c + (j-0.5)dt,
c      j=1,2,...,n, where dt = (d-c)/n.  n must be greater than 2.
c
c    nbdcnd
c      indicates the type of boundary conditions at theta = c
c      and theta = d.
c
c      = 0  if the solution is periodic in theta, i.e.
c           u(i,j) = u(i,n+j).
c
c      = 1  if the solution is specified at theta = c and theta = d
c           (see note below).
c
c      = 2  if the solution is specified at theta = c and the derivative
c           of the solution with respect to theta is specified at
c           theta = d (see note below).
c
c      = 3  if the derivative of the solution with respect to theta is
c           specified at theta = c and theta = d.
c
c      = 4  if the derivative of the solution with respect to theta is
c           specified at theta = c and the solution is specified at
c           theta = d (see note below).
c
c      note:  when nbdcnd = 1, 2, or 4, do not use mbdcnd = 5 or 6 (the
c      former indicates that the solution is specified at r =  0; the
c      latter indicates the solution is unspecified at r = 0).  use
c      instead mbdcnd = 1 or 2.
c
c    bdc
c      a one dimensional array of length m that specifies the boundary
c      values of the solution at theta = c.   when nbdcnd = 1 or 2,
c
c               bdc(i) = u(r(i),c) ,              i=1,2,...,m.
c
c      when nbdcnd = 3 or 4,
c
c               bdc(i) = (d/dtheta)u(r(i),c),     i=1,2,...,m.
c
c      when nbdcnd = 0, bdc is a dummy variable.
c
c    bdd
c      a one-dimensional array of length m that specifies the boundary
c      values of the solution at theta = d.  when nbdcnd = 1 or 4,
c
c               bdd(i) = u(r(i),d) ,              i=1,2,...,m.
c
c      when nbdcnd = 2 or 3,
c
c               bdd(i) = (d/dtheta)u(r(i),d) ,    i=1,2,...,m.
c
c      when nbdcnd = 0, bdd is a dummy variable.
c
c    elmbda
c      the constant lambda in the helmholtz equation.  if lambda is
c      greater than 0, a solution may not exist.  however, hstplr will
c      attempt to find a solution.
c
c    f
c      a two-dimensional array that specifies the values of the right
c      side of the helmholtz equation.  for i=1,2,...,m and j=1,2,...,n
c
c               f(i,j) = f(r(i),theta(j)) .
c
c      f must be dimensioned at least m x n.
c
c    idimf
c      the row (or first) dimension of the array f as it appears in the
c      program calling hstplr.  this parameter is used to specify the
c      variable dimension of f.  idimf must be at least m.
c
c    w
c      a one-dimensional array that must be provided by the user for
c      work space.  w may require up to 13m + 4n + m*int(log2(n))
c      locations.  the actual number of locations used is computed by
c      hstplr and is returned in the location w(1).
c
c
c             * * * * * *   on output   * * * * * *
c
c    f
c      contains the solution u(i,j) of the finite difference
c      approximation for the grid point (r(i),theta(j)) for
c      i=1,2,...,m, j=1,2,...,n.
c
c    pertrb
c      if a combination of periodic, derivative, or unspecified
c      boundary conditions is specified for a poisson equation
c      (lambda = 0), a solution may not exist.  pertrb is a con-
c      stant, calculated and subtracted from f, which ensures
c      that a solution exists.  hstplr then computes this
c      solution, which is a least squares solution to the
c      original approximation.  this solution plus any constant is also
c      a solution; hence, the solution is not unique.  the value of
c      pertrb should be small compared to the right side f.
c      otherwise, a solution is obtained to an essentially different
c      problem.  this comparison should always be made to insure that
c      a meaningful solution has been obtained.
c
c    ierror
c      an error flag that indicates invalid input parameters.
c      except for numbers 0 and 11, a solution is not attempted.
c
c      =  0  no error
c
c      =  1  a .lt. 0
c
c      =  2  a .ge. b
c
c      =  3  mbdcnd .lt. 1 or mbdcnd .gt. 6
c
c      =  4  c .ge. d
c
c      =  5  n .le. 2
c
c      =  6  nbdcnd .lt. 0 or nbdcnd .gt. 4
c
c      =  7  a = 0 and mbdcnd = 3 or 4
c
c      =  8  a .gt. 0 and mbdcnd .ge. 5
c
c      =  9  mbdcnd .ge. 5 and nbdcnd .ne. 0 or 3
c
c      = 10  idimf .lt. m
c
c      = 11  lambda .gt. 0
c
c      = 12  m .le. 2
c
c      since this is the only means of indicating a possibly
c      incorrect call to hstplr, the user should test ierror after
c      the call.
c
c    w
c      w(1) contains the required length of w.
c
c *long description:
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   bda(n),bdb(n),bdc(m),bdd(m),f(idimf,n),
c     arguments      w(see argument list)
c
c     latest         june 1, 1977
c     revision
c
c     subprograms    hstplr,poistg,postg2,genbun,poisd2,poisn2,poisp2,
c     required       cosgen,merge,trix,tri3,pimach
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
c     specialist     roland sweet
c
c     language       fortran
c
c     history        written by roland sweet at ncar in february, 1977
c
c     algorithm      this subroutine defines the finite-difference
c                    equations, incorporates boundary data, adjusts the
c                    right side when the system is singular and calls
c                    either poistg or genbun which solves the linear
c                    system of equations.
c
c     space          8265(decimal) = 20111(octal) locations on the
c     required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hstplr is roughly proportional
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
c                        32       1-6       1-4         56
c                        64       1-6       1-4        230
c
c     portability    american national standards institute fortran.
c                    the machine dependent constant pi is defined in
c                    function pimach.
c
c     required       cos
c     resident
c     routines
c
c     reference      schumann, u. and r. sweet,'a direct method for
c                    the solution of poisson's equation with neumann
c                    boundary conditions on a staggered grid of
c                    arbitrary size,' j. comp. phys. 20(1976),
c                    pp. 171-182.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c***references  u. schumann and r. sweet, a direct method for the
c                 solution of poisson's equation with neumann boundary
c                 conditions on a staggered grid of arbitrary size,
c                 journal of computational physics 20, (1976),
c                 pp. 171-182.
c***routines called  genbun, poistg
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  hstplr
c
c
      dimension       f(idimf,*)
      dimension       bda(*)     ,bdb(*)     ,bdc(*)     ,bdd(*)     ,
     1                w(*)
c***first executable statement  hstplr
      ierror = 0
      if (a .lt. 0.) ierror = 1
      if (a .ge. b) ierror = 2
      if (mbdcnd.le.0 .or. mbdcnd.ge.7) ierror = 3
      if (c .ge. d) ierror = 4
      if (n .le. 2) ierror = 5
      if (nbdcnd.lt.0 .or. nbdcnd.ge.5) ierror = 6
      if (a.eq.0. .and. (mbdcnd.eq.3 .or. mbdcnd.eq.4)) ierror = 7
      if (a.gt.0. .and. mbdcnd.ge.5) ierror = 8
      if (mbdcnd.ge.5 .and. nbdcnd.ne.0 .and. nbdcnd.ne.3) ierror = 9
      if (idimf .lt. m) ierror = 10
      if (m .le. 2) ierror = 12
      if (ierror .ne. 0) return
      deltar = (b-a)/m
      dlrsq = deltar**2
      deltht = (d-c)/n
      dlthsq = deltht**2
      np = nbdcnd+1
      isw = 1
      mb = mbdcnd
      if (a.eq.0. .and. mbdcnd.eq.2) mb = 6
c
c     define a,b,c coefficients in w-array.
c
      iwb = m
      iwc = iwb+m
      iwr = iwc+m
      do 101 i=1,m
         j = iwr+i
         w(j) = a+(i-0.5)*deltar
         w(i) = (a+(i-1)*deltar)/dlrsq
         k = iwc+i
         w(k) = (a+i*deltar)/dlrsq
         k = iwb+i
         w(k) = (elmbda-2./dlrsq)*w(j)
  101 continue
      do 103 i=1,m
         j = iwr+i
         a1 = w(j)
         do 102 j=1,n
            f(i,j) = a1*f(i,j)
  102    continue
  103 continue
c
c     enter boundary data for r-boundaries.
c
      go to (104,104,106,106,108,108),mb
  104 a1 = 2.*w(1)
      w(iwb+1) = w(iwb+1)-w(1)
      do 105 j=1,n
         f(1,j) = f(1,j)-a1*bda(j)
  105 continue
      go to 108
  106 a1 = deltar*w(1)
      w(iwb+1) = w(iwb+1)+w(1)
      do 107 j=1,n
         f(1,j) = f(1,j)+a1*bda(j)
  107 continue
  108 go to (109,111,111,109,109,111),mb
  109 a1 = 2.*w(iwr)
      w(iwc) = w(iwc)-w(iwr)
      do 110 j=1,n
         f(m,j) = f(m,j)-a1*bdb(j)
  110 continue
      go to 113
  111 a1 = deltar*w(iwr)
      w(iwc) = w(iwc)+w(iwr)
      do 112 j=1,n
         f(m,j) = f(m,j)-a1*bdb(j)
  112 continue
c
c     enter boundary data for theta-boundaries.
c
  113 a1 = 2./dlthsq
      go to (123,114,114,116,116),np
  114 do 115 i=1,m
         j = iwr+i
         f(i,1) = f(i,1)-a1*bdc(i)/w(j)
  115 continue
      go to 118
  116 a1 = 1./deltht
      do 117 i=1,m
         j = iwr+i
         f(i,1) = f(i,1)+a1*bdc(i)/w(j)
  117 continue
  118 a1 = 2./dlthsq
      go to (123,119,121,121,119),np
  119 do 120 i=1,m
         j = iwr+i
         f(i,n) = f(i,n)-a1*bdd(i)/w(j)
  120 continue
      go to 123
  121 a1 = 1./deltht
      do 122 i=1,m
         j = iwr+i
         f(i,n) = f(i,n)-a1*bdd(i)/w(j)
  122 continue
  123 continue
c
c     adjust right side of singular problems to insure existence of a
c     solution.
c
      pertrb = 0.
      if (elmbda) 133,125,124
  124 ierror = 11
      go to 133
  125 go to (133,133,126,133,133,126),mb
  126 go to (127,133,133,127,133),np
  127 continue
      isw = 2
      do 129 j=1,n
         do 128 i=1,m
            pertrb = pertrb+f(i,j)
  128    continue
  129 continue
      pertrb = pertrb/(m*n*0.5*(a+b))
      do 131 i=1,m
         j = iwr+i
         a1 = pertrb*w(j)
         do 130 j=1,n
            f(i,j) = f(i,j)-a1
  130    continue
  131 continue
      a2 = 0.
      do 132 j=1,n
         a2 = a2+f(1,j)
  132 continue
      a2 = a2/w(iwr+1)
  133 continue
c
c     multiply i-th equation through by  r(i)*deltht**2
c
      do 135 i=1,m
         j = iwr+i
         a1 = dlthsq*w(j)
         w(i) = a1*w(i)
         j = iwc+i
         w(j) = a1*w(j)
         j = iwb+i
         w(j) = a1*w(j)
         do 134 j=1,n
            f(i,j) = a1*f(i,j)
  134    continue
  135 continue
      lp = nbdcnd
      w(1) = 0.
      w(iwr) = 0.
c
c     call poistg or genbun to solve the system of equations.
c
      if (lp .eq. 0) go to 136
      call poistg (lp,n,1,m,w,w(iwb+1),w(iwc+1),idimf,f,ierr1,w(iwr+1))
      go to 137
  136 call genbun (lp,n,1,m,w,w(iwb+1),w(iwc+1),idimf,f,ierr1,w(iwr+1))
  137 continue
      w(1) = w(iwr+1)+3*m
      if (a.ne.0. .or. mbdcnd.ne.2 .or. isw.ne.2) go to 141
      a1 = 0.
      do 138 j=1,n
         a1 = a1+f(1,j)
  138 continue
      a1 = (a1-dlrsq*a2/16.)/n
      if (nbdcnd .eq. 3) a1 = a1+(bdd(1)-bdc(1))/(d-c)
      a1 = bda(1)-a1
      do 140 i=1,m
         do 139 j=1,n
            f(i,j) = f(i,j)+a1
  139    continue
  140 continue
  141 continue
      return
      end
