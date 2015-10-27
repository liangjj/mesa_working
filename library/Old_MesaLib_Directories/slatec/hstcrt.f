*deck hstcrt
      subroutine hstcrt (a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd,
     +   bdc, bdd, elmbda, f, idimf, pertrb, ierror, w)
c***begin prologue  hstcrt
c***purpose  solve the standard five-point finite difference
c            approximation on a staggered grid to the helmholtz equation
c            in cartesian coordinates.
c***library   slatec (fishpack)
c***category  i2b1a1a
c***type      single precision (hstcrt-s)
c***keywords  elliptic, fishpack, helmholtz, pde
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c      hstcrt solves the standard five-point finite difference
c      approximation on a staggered grid to the helmholtz equation in
c      cartesian coordinates
c
c      (d/dx)(du/dx) + (d/dy)(du/dy) + lambda*u = f(x,y)
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c    a,b
c      the range of x, i.e. a .le. x .le. b.  a must be less than b.
c
c    m
c      the number of grid points in the interval (a,b).  the grid points
c      in the x-direction are given by x(i) = a + (i-0.5)dx for
c      i=1,2,...,m where dx =(b-a)/m.  m must be greater than 2.
c
c    mbdcnd
c      indicates the type of boundary conditions at x = a and x = b.
c
c      = 0  if the solution is periodic in x,
c           u(m+i,j) = u(i,j).
c
c      = 1  if the solution is specified at x = a and x = b.
c
c      = 2  if the solution is specified at x = a and the derivative
c           of the solution with respect to x is specified at x = b.
c
c      = 3  if the derivative of the solution with respect to x is
c           specified at x = a  and x = b.
c
c      = 4  if the derivative of the solution with respect to x is
c           specified at x = a  and the solution is specified at x = b.
c
c    bda
c      a one-dimensional array of length n that specifies the boundary
c      values (if any) of the solution at x = a.  when mbdcnd = 1 or 2,
c
c               bda(j) = u(a,y(j)) ,          j=1,2,...,n.
c
c      when mbdcnd = 3 or 4,
c
c               bda(j) = (d/dx)u(a,y(j)) ,    j=1,2,...,n.
c
c    bdb
c      a one-dimensional array of length n that specifies the boundary
c      values of the solution at x = b.  when mbdcnd = 1 or 4
c
c               bdb(j) = u(b,y(j)) ,          j=1,2,...,n.
c
c      when mbdcnd = 2 or 3
c
c               bdb(j) = (d/dx)u(b,y(j)) ,    j=1,2,...,n.
c
c    c,d
c      the range of y, i.e. c .le. y .le. d.  c must be less
c      than d.
c
c    n
c      the number of unknowns in the interval (c,d).  the unknowns in
c      the y-direction are given by y(j) = c + (j-0.5)dy,
c      j=1,2,...,n, where dy = (d-c)/n.  n must be greater than 2.
c
c    nbdcnd
c      indicates the type of boundary conditions at y = c
c      and y = d.
c
c      = 0  if the solution is periodic in y, i.e.
c           u(i,j) = u(i,n+j).
c
c      = 1  if the solution is specified at y = c and y = d.
c
c      = 2  if the solution is specified at y = c and the derivative
c           of the solution with respect to y is specified at y = d.
c
c      = 3  if the derivative of the solution with respect to y is
c           specified at y = c and y = d.
c
c      = 4  if the derivative of the solution with respect to y is
c           specified at y = c and the solution is specified at y = d.
c
c    bdc
c      a one dimensional array of length m that specifies the boundary
c      values of the solution at y = c.   when nbdcnd = 1 or 2,
c
c               bdc(i) = u(x(i),c) ,              i=1,2,...,m.
c
c      when nbdcnd = 3 or 4,
c
c               bdc(i) = (d/dy)u(x(i),c),     i=1,2,...,m.
c
c      when nbdcnd = 0, bdc is a dummy variable.
c
c    bdd
c      a one-dimensional array of length m that specifies the boundary
c      values of the solution at y = d.  when nbdcnd = 1 or 4,
c
c               bdd(i) = u(x(i),d) ,              i=1,2,...,m.
c
c      when nbdcnd = 2 or 3,
c
c               bdd(i) = (d/dy)u(x(i),d) ,    i=1,2,...,m.
c
c      when nbdcnd = 0, bdd is a dummy variable.
c
c    elmbda
c      the constant lambda in the helmholtz equation.  if lambda is
c      greater than 0, a solution may not exist.  however, hstcrt will
c      attempt to find a solution.
c
c    f
c      a two-dimensional array that specifies the values of the right
c      side of the helmholtz equation.  for i=1,2,...,m and j=1,2,...,n
c
c               f(i,j) = f(x(i),y(j)) .
c
c      f must be dimensioned at least m x n.
c
c    idimf
c      the row (or first) dimension of the array f as it appears in the
c      program calling hstcrt.  this parameter is used to specify the
c      variable dimension of f.  idimf must be at least m.
c
c    w
c      a one-dimensional array that must be provided by the user for
c      work space.  w may require up to 13m + 4n + m*int(log2(n))
c      locations.  the actual number of locations used is computed by
c      hstcrt and is returned in the location w(1).
c
c
c             * * * * * *   on output   * * * * * *
c
c    f
c      contains the solution u(i,j) of the finite difference
c      approximation for the grid point (x(i),y(j)) for
c      i=1,2,...,m, j=1,2,...,n.
c
c    pertrb
c      if a combination of periodic or derivative boundary conditions is
c      specified for a poisson equation (lambda = 0), a solution may not
c      exist.  pertrb is a constant, calculated and subtracted from f,
c      which ensures that a solution exists.  hstcrt then computes this
c      solution, which is a least squares solution to the original
c      approximation.  this solution plus any constant is also a
c      solution; hence, the solution is not unique.  the value of pertrb
c      should be small compared to the right side f.  otherwise, a
c      solution is obtained to an essentially different problem.  this
c      comparison should always be made to insure that a meaningful
c      solution has been obtained.
c
c    ierror
c      an error flag that indicates invalid input parameters.
c       except for numbers 0 and  6, a solution is not attempted.
c
c      =  0  no error
c
c      =  1  a .ge. b
c
c      =  2  mbdcnd .lt. 0 or mbdcnd .gt. 4
c
c      =  3  c .ge. d
c
c      =  4  n .le. 2
c
c      =  5  nbdcnd .lt. 0 or nbdcnd .gt. 4
c
c      =  6  lambda .gt. 0
c
c      =  7  idimf .lt. m
c
c      =  8  m .le. 2
c
c      since this is the only means of indicating a possibly
c      incorrect call to hstcrt, the user should test ierror after
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
c     subprograms    hstcrt,poistg,postg2,genbun,poisd2,poisn2,poisp2,
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
c     history        written by roland sweet at ncar in january , 1977
c
c     algorithm      this subroutine defines the finite-difference
c                    equations, incorporates boundary data, adjusts the
c                    right side when the system is singular and calls
c                    either poistg or genbun which solves the linear
c                    system of equations.
c
c     space          8131(decimal) = 17703(octal) locations on the
c     required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hstcrt is roughly proportional
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
c                        32       1-4       1-4         56
c                        64       1-4       1-4        230
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
c***end prologue  hstcrt
c
c
      dimension       f(idimf,*) ,bda(*)     ,bdb(*)     ,bdc(*)     ,
     1                bdd(*)     ,w(*)
c***first executable statement  hstcrt
      ierror = 0
      if (a .ge. b) ierror = 1
      if (mbdcnd.lt.0 .or. mbdcnd.gt.4) ierror = 2
      if (c .ge. d) ierror = 3
      if (n .le. 2) ierror = 4
      if (nbdcnd.lt.0 .or. nbdcnd.gt.4) ierror = 5
      if (idimf .lt. m) ierror = 7
      if (m .le. 2) ierror = 8
      if (ierror .ne. 0) return
      nperod = nbdcnd
      mperod = 0
      if (mbdcnd .gt. 0) mperod = 1
      deltax = (b-a)/m
      twdelx = 1./deltax
      delxsq = 2./deltax**2
      deltay = (d-c)/n
      twdely = 1./deltay
      delysq = deltay**2
      twdysq = 2./delysq
      np = nbdcnd+1
      mp = mbdcnd+1
c
c     define the a,b,c coefficients in w-array.
c
      id2 = m
      id3 = id2+m
      id4 = id3+m
      s = (deltay/deltax)**2
      st2 = 2.*s
      do 101 i=1,m
         w(i) = s
         j = id2+i
         w(j) = -st2+elmbda*delysq
         j = id3+i
         w(j) = s
  101 continue
c
c     enter boundary data for x-boundaries.
c
      go to (111,102,102,104,104),mp
  102 do 103 j=1,n
         f(1,j) = f(1,j)-bda(j)*delxsq
  103 continue
      w(id2+1) = w(id2+1)-w(1)
      go to 106
  104 do 105 j=1,n
         f(1,j) = f(1,j)+bda(j)*twdelx
  105 continue
      w(id2+1) = w(id2+1)+w(1)
  106 go to (111,107,109,109,107),mp
  107 do 108 j=1,n
         f(m,j) = f(m,j)-bdb(j)*delxsq
  108 continue
      w(id3) = w(id3)-w(1)
      go to 111
  109 do 110 j=1,n
         f(m,j) = f(m,j)-bdb(j)*twdelx
  110 continue
      w(id3) = w(id3)+w(1)
  111 continue
c
c     enter boundary data for y-boundaries.
c
      go to (121,112,112,114,114),np
  112 do 113 i=1,m
         f(i,1) = f(i,1)-bdc(i)*twdysq
  113 continue
      go to 116
  114 do 115 i=1,m
         f(i,1) = f(i,1)+bdc(i)*twdely
  115 continue
  116 go to (121,117,119,119,117),np
  117 do 118 i=1,m
         f(i,n) = f(i,n)-bdd(i)*twdysq
  118 continue
      go to 121
  119 do 120 i=1,m
         f(i,n) = f(i,n)-bdd(i)*twdely
  120 continue
  121 continue
      do 123 i=1,m
         do 122 j=1,n
            f(i,j) = f(i,j)*delysq
  122    continue
  123 continue
      if (mperod .eq. 0) go to 124
      w(1) = 0.
      w(id4) = 0.
  124 continue
      pertrb = 0.
      if (elmbda) 133,126,125
  125 ierror = 6
      go to 133
  126 go to (127,133,133,127,133),mp
  127 go to (128,133,133,128,133),np
c
c     for singular problems must adjust data to insure that a solution
c     will exist.
c
  128 continue
      s = 0.
      do 130 j=1,n
         do 129 i=1,m
            s = s+f(i,j)
  129    continue
  130 continue
      pertrb = s/(m*n)
      do 132 j=1,n
         do 131 i=1,m
            f(i,j) = f(i,j)-pertrb
  131    continue
  132 continue
      pertrb = pertrb/delysq
c
c     solve the equation.
c
  133 continue
      if (nperod .eq. 0) go to 134
      call poistg (nperod,n,mperod,m,w(1),w(id2+1),w(id3+1),idimf,f,
     1             ierr1,w(id4+1))
      go to 135
  134 continue
      call genbun (nperod,n,mperod,m,w(1),w(id2+1),w(id3+1),idimf,f,
     1             ierr1,w(id4+1))
  135 continue
      w(1) = w(id4+1)+3*m
      return
      end
