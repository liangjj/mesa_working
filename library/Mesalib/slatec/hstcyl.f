*deck hstcyl
      subroutine hstcyl (a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd,
     +   bdc, bdd, elmbda, f, idimf, pertrb, ierror, w)
c***begin prologue  hstcyl
c***purpose  solve the standard five-point finite difference
c            approximation on a staggered grid to the modified
c            helmholtz equation in cylindrical coordinates.
c***library   slatec (fishpack)
c***category  i2b1a1a
c***type      single precision (hstcyl-s)
c***keywords  cylindrical, elliptic, fishpack, helmholtz, pde
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c      hstcyl solves the standard five-point finite difference
c      approximation on a staggered grid to the modified helmholtz
c      equation in cylindrical coordinates
c
c          (1/r)(d/dr)(r(du/dr)) + (d/dz)(du/dz)c
c                      + lambda*(1/r**2)*u = f(r,z)
c
c      this two-dimensional modified helmholtz equation results
c      from the fourier transform of a three-dimensional poisson
c      equation.
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
c      = 1  if the solution is specified at r = a (see note below) and
c           r = b.
c
c      = 2  if the solution is specified at r = a (see note below) and
c           the derivative of the solution with respect to r is
c           specified at r = b.
c
c      = 3  if the derivative of the solution with respect to r is
c           specified at r = a (see note below) and r = b.
c
c      = 4  if the derivative of the solution with respect to r is
c           specified at r = a (see note below) and the solution is
c           specified at r = b.
c
c      = 5  if the solution is unspecified at r = a = 0 and the solution
c           is specified at r = b.
c
c      = 6  if the solution is unspecified at r = a = 0 and the
c           derivative of the solution with respect to r is specified at
c           r = b.
c
c      note:  if a = 0, do not use mbdcnd = 1,2,3, or 4, but instead
c             use mbdcnd = 5 or 6.  the resulting approximation gives
c             the only meaningful boundary condition, i.e. du/dr = 0.
c             (see d. greenspan, 'introductory numerical analysis of
c             elliptic boundary value problems,' harper and row, 1965,
c             chapter 5.)
c
c    bda
c      a one-dimensional array of length n that specifies the boundary
c      values (if any) of the solution at r = a.  when mbdcnd = 1 or 2,
c
c               bda(j) = u(a,z(j)) ,          j=1,2,...,n.
c
c      when mbdcnd = 3 or 4,
c
c               bda(j) = (d/dr)u(a,z(j)) ,    j=1,2,...,n.
c
c      when mbdcnd = 5 or 6, bda is a dummy variable.
c
c    bdb
c      a one-dimensional array of length n that specifies the boundary
c      values of the solution at r = b.  when mbdcnd = 1,4, or 5,
c
c               bdb(j) = u(b,z(j)) ,          j=1,2,...,n.
c
c      when mbdcnd = 2,3, or 6,
c
c               bdb(j) = (d/dr)u(b,z(j)) ,    j=1,2,...,n.
c
c    c,d
c      the range of z, i.e. c .le. z .le. d.  c must be less
c      than d.
c
c    n
c      the number of unknowns in the interval (c,d).  the unknowns in
c      the z-direction are given by z(j) = c + (j-0.5)dz,
c      j=1,2,...,n, where dz = (d-c)/n.  n must be greater than 2.
c
c    nbdcnd
c      indicates the type of boundary conditions at z = c
c      and z = d.
c
c      = 0  if the solution is periodic in z, i.e.
c           u(i,j) = u(i,n+j).
c
c      = 1  if the solution is specified at z = c and z = d.
c
c      = 2  if the solution is specified at z = c and the derivative
c           of the solution with respect to z is specified at
c           z = d.
c
c      = 3  if the derivative of the solution with respect to z is
c           specified at z = c and z = d.
c
c      = 4  if the derivative of the solution with respect to z is
c           specified at z = c and the solution is specified at
c           z = d.
c
c    bdc
c      a one dimensional array of length m that specifies the boundary
c      values of the solution at z = c.   when nbdcnd = 1 or 2,
c
c               bdc(i) = u(r(i),c) ,              i=1,2,...,m.
c
c      when nbdcnd = 3 or 4,
c
c               bdc(i) = (d/dz)u(r(i),c),         i=1,2,...,m.
c
c      when nbdcnd = 0, bdc is a dummy variable.
c
c    bdd
c      a one-dimensional array of length m that specifies the boundary
c      values of the solution at z = d.  when nbdcnd = 1 or 4,
c
c               bdd(i) = u(r(i),d) ,              i=1,2,...,m.
c
c      when nbdcnd = 2 or 3,
c
c               bdd(i) = (d/dz)u(r(i),d) ,        i=1,2,...,m.
c
c      when nbdcnd = 0, bdd is a dummy variable.
c
c    elmbda
c      the constant lambda in the modified helmholtz equation.  if
c      lambda is greater than 0, a solution may not exist.  however,
c      hstcyl will attempt to find a solution.  lambda must be zero
c      when mbdcnd = 5 or 6.
c
c    f
c      a two-dimensional array that specifies the values of the right
c      side of the modified helmholtz equation.  for i=1,2,...,m
c      and j=1,2,...,n
c
c               f(i,j) = f(r(i),z(j)) .
c
c      f must be dimensioned at least m x n.
c
c    idimf
c      the row (or first) dimension of the array f as it appears in the
c      program calling hstcyl.  this parameter is used to specify the
c      variable dimension of f.  idimf must be at least m.
c
c    w
c      a one-dimensional array that must be provided by the user for
c      work space.  w may require up to 13m + 4n + m*int(log2(n))
c      locations.  the actual number of locations used is computed by
c      hstcyl and is returned in the location w(1).
c
c
c             * * * * * *   on output   * * * * * *
c
c    f
c      contains the solution u(i,j) of the finite difference
c      approximation for the grid point (r(i),z(j)) for
c      i=1,2,...,m, j=1,2,...,n.
c
c    pertrb
c      if a combination of periodic, derivative, or unspecified
c      boundary conditions is specified for a poisson equation
c      (lambda = 0), a solution may not exist.  pertrb is a con-
c      stant, calculated and subtracted from f, which ensures
c      that a solution exists.  hstcyl then computes this
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
c      =  7  a = 0 and mbdcnd = 1,2,3, or 4
c
c      =  8  a .gt. 0 and mbdcnd .ge. 5
c
c      =  9  m .le. 2
c
c      = 10  idimf .lt. m
c
c      = 11  lambda .gt. 0
c
c      = 12  a=0, mbdcnd .ge. 5, elmbda .ne. 0
c
c      since this is the only means of indicating a possibly
c      incorrect call to hstcyl, the user should test ierror after
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
c     subprograms    hstcyl,poistg,postg2,genbun,poisd2,poisn2,poisp2,
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
c     history        written by roland sweet at ncar in march, 1977
c
c     algorithm      this subroutine defines the finite-difference
c                    equations, incorporates boundary data, adjusts the
c                    right side when the system is singular and calls
c                    either poistg or genbun which solves the linear
c                    system of equations.
c
c     space          8228(decimal) = 20044(octal) locations on the
c     required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hstcyl is roughly proportional
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
c***end prologue  hstcyl
c
c
      dimension       f(idimf,*) ,bda(*)     ,bdb(*)     ,bdc(*)     ,
     1                bdd(*)     ,w(*)
c***first executable statement  hstcyl
      ierror = 0
      if (a .lt. 0.) ierror = 1
      if (a .ge. b) ierror = 2
      if (mbdcnd.le.0 .or. mbdcnd.ge.7) ierror = 3
      if (c .ge. d) ierror = 4
      if (n .le. 2) ierror = 5
      if (nbdcnd.lt.0 .or. nbdcnd.ge.5) ierror = 6
      if (a.eq.0. .and. mbdcnd.ne.5 .and. mbdcnd.ne.6) ierror = 7
      if (a.gt.0. .and. mbdcnd.ge.5) ierror = 8
      if (idimf .lt. m) ierror = 10
      if (m .le. 2) ierror = 9
      if (a.eq.0. .and. mbdcnd.ge.5 .and. elmbda.ne.0.) ierror = 12
      if (ierror .ne. 0) return
      deltar = (b-a)/m
      dlrsq = deltar**2
      deltht = (d-c)/n
      dlthsq = deltht**2
      np = nbdcnd+1
c
c     define a,b,c coefficients in w-array.
c
      iwb = m
      iwc = iwb+m
      iwr = iwc+m
      do 101 i=1,m
         j = iwr+i
         w(j) = a+(i-0.5)*deltar
         w(i) = (a+(i-1)*deltar)/(dlrsq*w(j))
         k = iwc+i
         w(k) = (a+i*deltar)/(dlrsq*w(j))
         k = iwb+i
         w(k) = elmbda/w(j)**2-2./dlrsq
  101 continue
c
c     enter boundary data for r-boundaries.
c
      go to (102,102,104,104,106,106),mbdcnd
  102 a1 = 2.*w(1)
      w(iwb+1) = w(iwb+1)-w(1)
      do 103 j=1,n
         f(1,j) = f(1,j)-a1*bda(j)
  103 continue
      go to 106
  104 a1 = deltar*w(1)
      w(iwb+1) = w(iwb+1)+w(1)
      do 105 j=1,n
         f(1,j) = f(1,j)+a1*bda(j)
  105 continue
  106 continue
      go to (107,109,109,107,107,109),mbdcnd
  107 w(iwc) = w(iwc)-w(iwr)
      a1 = 2.*w(iwr)
      do 108 j=1,n
         f(m,j) = f(m,j)-a1*bdb(j)
  108 continue
      go to 111
  109 w(iwc) = w(iwc)+w(iwr)
      a1 = deltar*w(iwr)
      do 110 j=1,n
         f(m,j) = f(m,j)-a1*bdb(j)
  110 continue
c
c     enter boundary data for theta-boundaries.
c
  111 a1 = 2./dlthsq
      go to (121,112,112,114,114),np
  112 do 113 i=1,m
         f(i,1) = f(i,1)-a1*bdc(i)
  113 continue
      go to 116
  114 a1 = 1./deltht
      do 115 i=1,m
         f(i,1) = f(i,1)+a1*bdc(i)
  115 continue
  116 a1 = 2./dlthsq
      go to (121,117,119,119,117),np
  117 do 118 i=1,m
         f(i,n) = f(i,n)-a1*bdd(i)
  118 continue
      go to 121
  119 a1 = 1./deltht
      do 120 i=1,m
         f(i,n) = f(i,n)-a1*bdd(i)
  120 continue
  121 continue
c
c     adjust right side of singular problems to insure existence of a
c     solution.
c
      pertrb = 0.
      if (elmbda) 130,123,122
  122 ierror = 11
      go to 130
  123 go to (130,130,124,130,130,124),mbdcnd
  124 go to (125,130,130,125,130),np
  125 continue
      do 127 i=1,m
         a1 = 0.
         do 126 j=1,n
            a1 = a1+f(i,j)
  126    continue
         j = iwr+i
         pertrb = pertrb+a1*w(j)
  127 continue
      pertrb = pertrb/(m*n*0.5*(a+b))
      do 129 i=1,m
         do 128 j=1,n
            f(i,j) = f(i,j)-pertrb
  128    continue
  129 continue
  130 continue
c
c     multiply i-th equation through by  deltht**2
c
      do 132 i=1,m
         w(i) = w(i)*dlthsq
         j = iwc+i
         w(j) = w(j)*dlthsq
         j = iwb+i
         w(j) = w(j)*dlthsq
         do 131 j=1,n
            f(i,j) = f(i,j)*dlthsq
  131    continue
  132 continue
      lp = nbdcnd
      w(1) = 0.
      w(iwr) = 0.
c
c     call genbun to solve the system of equations.
c
      if (nbdcnd .eq. 0) go to 133
      call poistg (lp,n,1,m,w,w(iwb+1),w(iwc+1),idimf,f,ierr1,w(iwr+1))
      go to 134
  133 call genbun (lp,n,1,m,w,w(iwb+1),w(iwc+1),idimf,f,ierr1,w(iwr+1))
  134 continue
      w(1) = w(iwr+1)+3*m
      return
      end
