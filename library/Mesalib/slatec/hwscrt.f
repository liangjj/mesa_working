*deck hwscrt
      subroutine hwscrt (a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd,
     +   bdc, bdd, elmbda, f, idimf, pertrb, ierror, w)
c***begin prologue  hwscrt
c***purpose  solves the standard five-point finite difference
c            approximation to the helmholtz equation in cartesian
c            coordinates.
c***library   slatec (fishpack)
c***category  i2b1a1a
c***type      single precision (hwscrt-s)
c***keywords  cartesian, elliptic, fishpack, helmholtz, pde
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c     subroutine hwscrt solves the standard five-point finite
c     difference approximation to the helmholtz equation in cartesian
c     coordinates:
c
c          (d/dx)(du/dx) + (d/dy)(du/dy) + lambda*u = f(x,y).
c
c
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c     a,b
c       the range of x, i.e., a .le. x .le. b.  a must be less than b.
c
c     m
c       the number of panels into which the interval (a,b) is
c       subdivided.  hence, there will be m+1 grid points in the
c       x-direction given by x(i) = a+(i-1)dx for i = 1,2,...,m+1,
c       where dx = (b-a)/m is the panel width. m must be greater than 3.
c
c     mbdcnd
c       indicates the type of boundary conditions at x = a and x = b.
c
c       = 0  if the solution is periodic in x, i.e., u(i,j) = u(m+i,j).
c       = 1  if the solution is specified at x = a and x = b.
c       = 2  if the solution is specified at x = a and the derivative of
c            the solution with respect to x is specified at x = b.
c       = 3  if the derivative of the solution with respect to x is
c            specified at x = a and x = b.
c       = 4  if the derivative of the solution with respect to x is
c            specified at x = a and the solution is specified at x = b.
c
c     bda
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to x at x = a.
c       when mbdcnd = 3 or 4,
c
c            bda(j) = (d/dx)u(a,y(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bda is a dummy variable.
c
c     bdb
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to x at x = b.
c       when mbdcnd = 2 or 3,
c
c            bdb(j) = (d/dx)u(b,y(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value bdb is a dummy variable.
c
c     c,d
c       the range of y, i.e., c .le. y .le. d.  c must be less than d.
c
c     n
c       the number of panels into which the interval (c,d) is
c       subdivided.  hence, there will be n+1 grid points in the
c       y-direction given by y(j) = c+(j-1)dy for j = 1,2,...,n+1, where
c       dy = (d-c)/n is the panel width.  n must be greater than 3.
c
c     nbdcnd
c       indicates the type of boundary conditions at y = c and y = d.
c
c       = 0  if the solution is periodic in y, i.e., u(i,j) = u(i,n+j).
c       = 1  if the solution is specified at y = c and y = d.
c       = 2  if the solution is specified at y = c and the derivative of
c            the solution with respect to y is specified at y = d.
c       = 3  if the derivative of the solution with respect to y is
c            specified at y = c and y = d.
c       = 4  if the derivative of the solution with respect to y is
c            specified at y = c and the solution is specified at y = d.
c
c     bdc
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to y at y = c.
c       when nbdcnd = 3 or 4,
c
c            bdc(i) = (d/dy)u(x(i),c), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdc is a dummy variable.
c
c     bdd
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to y at y = d.
c       when nbdcnd = 2 or 3,
c
c            bdd(i) = (d/dy)u(x(i),d), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdd is a dummy variable.
c
c     elmbda
c       the constant lambda in the helmholtz equation.  if
c       lambda .gt. 0, a solution may not exist.  however, hwscrt will
c       attempt to find a solution.
c
c     f
c       a two-dimensional array which specifies the values of the right
c       side of the helmholtz equation and boundary values (if any).
c       for i = 2,3,...,m and j = 2,3,...,n
c
c            f(i,j) = f(x(i),y(j)).
c
c       on the boundaries f is defined by
c
c            mbdcnd     f(1,j)        f(m+1,j)
c            ------     ---------     --------
c
c              0        f(a,y(j))     f(a,y(j))
c              1        u(a,y(j))     u(b,y(j))
c              2        u(a,y(j))     f(b,y(j))     j = 1,2,...,n+1
c              3        f(a,y(j))     f(b,y(j))
c              4        f(a,y(j))     u(b,y(j))
c
c
c            nbdcnd     f(i,1)        f(i,n+1)
c            ------     ---------     --------
c
c              0        f(x(i),c)     f(x(i),c)
c              1        u(x(i),c)     u(x(i),d)
c              2        u(x(i),c)     f(x(i),d)     i = 1,2,...,m+1
c              3        f(x(i),c)     f(x(i),d)
c              4        f(x(i),c)     u(x(i),d)
c
c       f must be dimensioned at least (m+1)*(n+1).
c
c       note:
c
c       if the table calls for both the solution u and the right side f
c       at a corner then the solution must be specified.
c
c     idimf
c       the row (or first) dimension of the array f as it appears in the
c       program calling hwscrt.  this parameter is used to specify the
c       variable dimension of f.  idimf must be at least m+1  .
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space.  w may require up to 4*(n+1) +
c       (13 + int(log2(n+1)))*(m+1) locations.  the actual number of
c       locations used is computed by hwscrt and is returned in location
c       w(1).
c
c
c             * * * * * *   on output     * * * * * *
c
c     f
c       contains the solution u(i,j) of the finite difference
c       approximation for the grid point (x(i),y(j)), i = 1,2,...,m+1,
c       j = 1,2,...,n+1  .
c
c     pertrb
c       if a combination of periodic or derivative boundary conditions
c       is specified for a poisson equation (lambda = 0), a solution may
c       not exist.  pertrb is a constant, calculated and subtracted from
c       f, which ensures that a solution exists.  hwscrt then computes
c       this solution, which is a least squares solution to the original
c       approximation.  this solution plus any constant is also a
c       solution.  hence, the solution is not unique.  the value of
c       pertrb should be small compared to the right side f.  otherwise,
c       a solution is obtained to an essentially different problem.
c       this comparison should always be made to insure that a
c       meaningful solution has been obtained.
c
c     ierror
c       an error flag that indicates invalid input parameters.  except
c       for numbers 0 and 6, a solution is not attempted.
c
c       = 0  no error.
c       = 1  a .ge. b.
c       = 2  mbdcnd .lt. 0 or mbdcnd .gt. 4  .
c       = 3  c .ge. d.
c       = 4  n .le. 3
c       = 5  nbdcnd .lt. 0 or nbdcnd .gt. 4  .
c       = 6  lambda .gt. 0  .
c       = 7  idimf .lt. m+1  .
c       = 8  m .le. 3
c
c       since this is the only means of indicating a possibly incorrect
c       call to hwscrt, the user should test ierror after the call.
c
c     w
c       w(1) contains the required length of w.
c
c *long description:
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c
c     dimension of   bda(n+1),bdb(n+1),bdc(m+1),bdd(m+1),f(idimf,n+1),
c     arguments      w(see argument list)
c
c     latest         june 1, 1976
c     revision
c
c     subprograms    hwscrt,genbun,poisd2,poisn2,poisp2,cosgen,merge,
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
c     specialist     roland sweet
c
c     language       fortran
c
c     history        standardized september 1, 1973
c                    revised april 1, 1976
c
c     algorithm      the routine defines the finite difference
c                    equations, incorporates boundary data, and adjusts
c                    the right side of singular systems and then calls
c                    genbun to solve the system.
c
c     space          13110(octal) = 5704(decimal) locations on the ncar
c     required       control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hwscrt is roughly proportional
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
c     reference      swarztrauber, p. and r. sweet, 'efficient fortran
c                    subprograms for the solution of elliptic equations'
c                    ncar tn/ia-109, july, 1975, 138 pp.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c***references  p. n. swarztrauber and r. sweet, efficient fortran
c                 subprograms for the solution of elliptic equations,
c                 ncar tn/ia-109, july 1975, 138 pp.
c***routines called  genbun
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  hwscrt
c
c
      dimension       f(idimf,*)
      dimension       bda(*)     ,bdb(*)     ,bdc(*)     ,bdd(*)     ,
     1                w(*)
c***first executable statement  hwscrt
      ierror = 0
      if (a .ge. b) ierror = 1
      if (mbdcnd.lt.0 .or. mbdcnd.gt.4) ierror = 2
      if (c .ge. d) ierror = 3
      if (n .le. 3) ierror = 4
      if (nbdcnd.lt.0 .or. nbdcnd.gt.4) ierror = 5
      if (idimf .lt. m+1) ierror = 7
      if (m .le. 3) ierror = 8
      if (ierror .ne. 0) return
      nperod = nbdcnd
      mperod = 0
      if (mbdcnd .gt. 0) mperod = 1
      deltax = (b-a)/m
      twdelx = 2./deltax
      delxsq = 1./deltax**2
      deltay = (d-c)/n
      twdely = 2./deltay
      delysq = 1./deltay**2
      np = nbdcnd+1
      np1 = n+1
      mp = mbdcnd+1
      mp1 = m+1
      nstart = 1
      nstop = n
      nskip = 1
      go to (104,101,102,103,104),np
  101 nstart = 2
      go to 104
  102 nstart = 2
  103 nstop = np1
      nskip = 2
  104 nunk = nstop-nstart+1
c
c     enter boundary data for x-boundaries.
c
      mstart = 1
      mstop = m
      mskip = 1
      go to (117,105,106,109,110),mp
  105 mstart = 2
      go to 107
  106 mstart = 2
      mstop = mp1
      mskip = 2
  107 do 108 j=nstart,nstop
         f(2,j) = f(2,j)-f(1,j)*delxsq
  108 continue
      go to 112
  109 mstop = mp1
      mskip = 2
  110 do 111 j=nstart,nstop
         f(1,j) = f(1,j)+bda(j)*twdelx
  111 continue
  112 go to (113,115),mskip
  113 do 114 j=nstart,nstop
         f(m,j) = f(m,j)-f(mp1,j)*delxsq
  114 continue
      go to 117
  115 do 116 j=nstart,nstop
         f(mp1,j) = f(mp1,j)-bdb(j)*twdelx
  116 continue
  117 munk = mstop-mstart+1
c
c     enter boundary data for y-boundaries.
c
      go to (127,118,118,120,120),np
  118 do 119 i=mstart,mstop
         f(i,2) = f(i,2)-f(i,1)*delysq
  119 continue
      go to 122
  120 do 121 i=mstart,mstop
         f(i,1) = f(i,1)+bdc(i)*twdely
  121 continue
  122 go to (123,125),nskip
  123 do 124 i=mstart,mstop
         f(i,n) = f(i,n)-f(i,np1)*delysq
  124 continue
      go to 127
  125 do 126 i=mstart,mstop
         f(i,np1) = f(i,np1)-bdd(i)*twdely
  126 continue
c
c    multiply right side by deltay**2.
c
  127 delysq = deltay*deltay
      do 129 i=mstart,mstop
         do 128 j=nstart,nstop
            f(i,j) = f(i,j)*delysq
  128    continue
  129 continue
c
c     define the a,b,c coefficients in w-array.
c
      id2 = munk
      id3 = id2+munk
      id4 = id3+munk
      s = delysq*delxsq
      st2 = 2.*s
      do 130 i=1,munk
         w(i) = s
         j = id2+i
         w(j) = -st2+elmbda*delysq
         j = id3+i
         w(j) = s
  130 continue
      if (mp .eq. 1) go to 131
      w(1) = 0.
      w(id4) = 0.
  131 continue
      go to (135,135,132,133,134),mp
  132 w(id2) = st2
      go to 135
  133 w(id2) = st2
  134 w(id3+1) = st2
  135 continue
      pertrb = 0.
      if (elmbda) 144,137,136
  136 ierror = 6
      go to 144
  137 if ((nbdcnd.eq.0 .or. nbdcnd.eq.3) .and.
     1    (mbdcnd.eq.0 .or. mbdcnd.eq.3)) go to 138
      go to 144
c
c     for singular problems must adjust data to insure that a solution
c     will exist.
c
  138 a1 = 1.
      a2 = 1.
      if (nbdcnd .eq. 3) a2 = 2.
      if (mbdcnd .eq. 3) a1 = 2.
      s1 = 0.
      msp1 = mstart+1
      mstm1 = mstop-1
      nsp1 = nstart+1
      nstm1 = nstop-1
      do 140 j=nsp1,nstm1
         s = 0.
         do 139 i=msp1,mstm1
            s = s+f(i,j)
  139    continue
         s1 = s1+s*a1+f(mstart,j)+f(mstop,j)
  140 continue
      s1 = a2*s1
      s = 0.
      do 141 i=msp1,mstm1
         s = s+f(i,nstart)+f(i,nstop)
  141 continue
      s1 = s1+s*a1+f(mstart,nstart)+f(mstart,nstop)+f(mstop,nstart)+
     1     f(mstop,nstop)
      s = (2.+(nunk-2)*a2)*(2.+(munk-2)*a1)
      pertrb = s1/s
      do 143 j=nstart,nstop
         do 142 i=mstart,mstop
            f(i,j) = f(i,j)-pertrb
  142    continue
  143 continue
      pertrb = pertrb/delysq
c
c     solve the equation.
c
  144 call genbun (nperod,nunk,mperod,munk,w(1),w(id2+1),w(id3+1),
     1             idimf,f(mstart,nstart),ierr1,w(id4+1))
      w(1) = w(id4+1)+3*munk
c
c     fill in identical values when have periodic boundary conditions.
c
      if (nbdcnd .ne. 0) go to 146
      do 145 i=mstart,mstop
         f(i,np1) = f(i,1)
  145 continue
  146 if (mbdcnd .ne. 0) go to 148
      do 147 j=nstart,nstop
         f(mp1,j) = f(1,j)
  147 continue
      if (nbdcnd .eq. 0) f(mp1,np1) = f(1,np1)
  148 continue
      return
      end
