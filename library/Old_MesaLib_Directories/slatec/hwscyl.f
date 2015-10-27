*deck hwscyl
      subroutine hwscyl (a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd,
     +   bdc, bdd, elmbda, f, idimf, pertrb, ierror, w)
c***begin prologue  hwscyl
c***purpose  solve a standard finite difference approximation
c            to the helmholtz equation in cylindrical coordinates.
c***library   slatec (fishpack)
c***category  i2b1a1a
c***type      single precision (hwscyl-s)
c***keywords  cylindrical, elliptic, fishpack, helmholtz, pde
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c     subroutine hwscyl solves a finite difference approximation to the
c     helmholtz equation in cylindrical coordinates:
c
c          (1/r)(d/dr)(r(du/dr)) + (d/dz)(du/dz)
c
c                                + (lambda/r**2)u = f(r,z)
c
c     this modified helmholtz equation results from the fourier
c     transform of the three-dimensional poisson equation.
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c     a,b
c       the range of r, i.e., a .le. r .le. b.  a must be less than b
c       and a must be non-negative.
c
c     m
c       the number of panels into which the interval (a,b) is
c       subdivided.  hence, there will be m+1 grid points in the
c       r-direction given by r(i) = a+(i-1)dr, for i = 1,2,...,m+1,
c       where dr = (b-a)/m is the panel width. m must be greater than 3.
c
c     mbdcnd
c       indicates the type of boundary conditions at r = a and r = b.
c
c       = 1  if the solution is specified at r = a and r = b.
c       = 2  if the solution is specified at r = a and the derivative of
c            the solution with respect to r is specified at r = b.
c       = 3  if the derivative of the solution with respect to r is
c            specified at r = a (see note below) and r = b.
c       = 4  if the derivative of the solution with respect to r is
c            specified at r = a (see note below) and the solution is
c            specified at r = b.
c       = 5  if the solution is unspecified at r = a = 0 and the
c            solution is specified at r = b.
c       = 6  if the solution is unspecified at r = a = 0 and the
c            derivative of the solution with respect to r is specified
c            at r = b.
c
c       note:  if a = 0, do not use mbdcnd = 3 or 4, but instead use
c              mbdcnd = 1,2,5, or 6  .
c
c     bda
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to r at r = a.
c       when mbdcnd = 3 or 4,
c
c            bda(j) = (d/dr)u(a,z(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bda is a dummy variable.
c
c     bdb
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to r at r = b.
c       when mbdcnd = 2,3, or 6,
c
c            bdb(j) = (d/dr)u(b,z(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bdb is a dummy variable.
c
c     c,d
c       the range of z, i.e., c .le. z .le. d.  c must be less than d.
c
c     n
c       the number of panels into which the interval (c,d) is
c       subdivided.  hence, there will be n+1 grid points in the
c       z-direction given by z(j) = c+(j-1)dz, for j = 1,2,...,n+1,
c       where dz = (d-c)/n is the panel width. n must be greater than 3.
c
c     nbdcnd
c       indicates the type of boundary conditions at z = c and z = d.
c
c       = 0  if the solution is periodic in z, i.e., u(i,1) = u(i,n+1).
c       = 1  if the solution is specified at z = c and z = d.
c       = 2  if the solution is specified at z = c and the derivative of
c            the solution with respect to z is specified at z = d.
c       = 3  if the derivative of the solution with respect to z is
c            specified at z = c and z = d.
c       = 4  if the derivative of the solution with respect to z is
c            specified at z = c and the solution is specified at z = d.
c
c     bdc
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to z at z = c.
c       when nbdcnd = 3 or 4,
c
c            bdc(i) = (d/dz)u(r(i),c), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdc is a dummy variable.
c
c     bdd
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to z at z = d.
c       when nbdcnd = 2 or 3,
c
c            bdd(i) = (d/dz)u(r(i),d), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdd is a dummy variable.
c
c     elmbda
c       the constant lambda in the helmholtz equation.  if
c       lambda .gt. 0, a solution may not exist.  however, hwscyl will
c       attempt to find a solution.  lambda must be zero when
c       mbdcnd = 5 or 6  .
c
c     f
c       a two-dimensional array that specifies the values of the right
c       side of the helmholtz equation and boundary data (if any).  for
c       i = 2,3,...,m and j = 2,3,...,n
c
c            f(i,j) = f(r(i),z(j)).
c
c       on the boundaries f is defined by
c
c            mbdcnd   f(1,j)            f(m+1,j)
c            ------   ---------         ---------
c
c              1      u(a,z(j))         u(b,z(j))
c              2      u(a,z(j))         f(b,z(j))
c              3      f(a,z(j))         f(b,z(j))   j = 1,2,...,n+1
c              4      f(a,z(j))         u(b,z(j))
c              5      f(0,z(j))         u(b,z(j))
c              6      f(0,z(j))         f(b,z(j))
c
c            nbdcnd   f(i,1)            f(i,n+1)
c            ------   ---------         ---------
c
c              0      f(r(i),c)         f(r(i),c)
c              1      u(r(i),c)         u(r(i),d)
c              2      u(r(i),c)         f(r(i),d)   i = 1,2,...,m+1
c              3      f(r(i),c)         f(r(i),d)
c              4      f(r(i),c)         u(r(i),d)
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
c       program calling hwscyl.  this parameter is used to specify the
c       variable dimension of f.  idimf must be at least m+1  .
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space.  w may require up to 4*(n+1) +
c       (13 + int(log2(n+1)))*(m+1) locations.  the actual number of
c       locations used is computed by hwscyl and is returned in location
c       w(1).
c
c
c             * * * * * *   on output     * * * * * *
c
c     f
c       contains the solution u(i,j) of the finite difference
c       approximation for the grid point (r(i),z(j)), i = 1,2,...,m+1,
c       j = 1,2,...,n+1  .
c
c     pertrb
c       if one specifies a combination of periodic, derivative, and
c       unspecified boundary conditions for a poisson equation
c       (lambda = 0), a solution may not exist.  pertrb is a constant,
c       calculated and subtracted from f, which ensures that a solution
c       exists.  hwscyl then computes this solution, which is a least
c       squares solution to the original approximation.  this solution
c       plus any constant is also a solution.  hence, the solution is
c       not unique.  the value of pertrb should be small compared to the
c       right side f.  otherwise, a solution is obtained to an
c       essentially different problem.  this comparison should always
c       be made to insure that a meaningful solution has been obtained.
c
c     ierror
c       an error flag which indicates invalid input parameters.  except
c       for numbers 0 and 11, a solution is not attempted.
c
c       =  0  no error.
c       =  1  a .lt. 0  .
c       =  2  a .ge. b.
c       =  3  mbdcnd .lt. 1 or mbdcnd .gt. 6  .
c       =  4  c .ge. d.
c       =  5  n .le. 3
c       =  6  nbdcnd .lt. 0 or nbdcnd .gt. 4  .
c       =  7  a = 0, mbdcnd = 3 or 4  .
c       =  8  a .gt. 0, mbdcnd .ge. 5  .
c       =  9  a = 0, lambda .ne. 0, mbdcnd .ge. 5  .
c       = 10  idimf .lt. m+1  .
c       = 11  lambda .gt. 0  .
c       = 12  m .le. 3
c
c       since this is the only means of indicating a possibly incorrect
c       call to hwscyl, the user should test ierror after the call.
c
c     w
c       w(1) contains the required length of w.
c
c *long description:
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   bda(n+1),bdb(n+1),bdc(m+1),bdd(m+1),f(idimf,n+1),
c     arguments      w(see argument list)
c
c     latest         june 1, 1976
c     revision
c
c     subprograms    hwscyl,genbun,poisd2,poisn2,poisp2,cosgen,merge,
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
c     space          5818(decimal) = 13272(octal) locations on the ncar
c     required       control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hwscyl is roughly proportional
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
c                        32        1         0          31
c                        32        1         1          23
c                        32        3         3          36
c                        64        1         0         128
c                        64        1         1          96
c                        64        3         3         142
c
c     portability    american national standards institute fortran.
c                    the machine dependent constant pi is defined in
c                    function pimach.
c
c     required       cos
c     resident
c     routines
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
c***end prologue  hwscyl
c
c
      dimension       f(idimf,*)
      dimension       bda(*)     ,bdb(*)     ,bdc(*)     ,bdd(*)     ,
     1                w(*)
c***first executable statement  hwscyl
      ierror = 0
      if (a .lt. 0.) ierror = 1
      if (a .ge. b) ierror = 2
      if (mbdcnd.le.0 .or. mbdcnd.ge.7) ierror = 3
      if (c .ge. d) ierror = 4
      if (n .le. 3) ierror = 5
      if (nbdcnd.le.-1 .or. nbdcnd.ge.5) ierror = 6
      if (a.eq.0. .and. (mbdcnd.eq.3 .or. mbdcnd.eq.4)) ierror = 7
      if (a.gt.0. .and. mbdcnd.ge.5) ierror = 8
      if (a.eq.0. .and. elmbda.ne.0. .and. mbdcnd.ge.5) ierror = 9
      if (idimf .lt. m+1) ierror = 10
      if (m .le. 3) ierror = 12
      if (ierror .ne. 0) return
      mp1 = m+1
      deltar = (b-a)/m
      dlrby2 = deltar/2.
      dlrsq = deltar**2
      np1 = n+1
      deltht = (d-c)/n
      dlthsq = deltht**2
      np = nbdcnd+1
c
c     define range of indices i and j for unknowns u(i,j).
c
      mstart = 2
      mstop = m
      go to (104,103,102,101,101,102),mbdcnd
  101 mstart = 1
      go to 104
  102 mstart = 1
  103 mstop = mp1
  104 munk = mstop-mstart+1
      nstart = 1
      nstop = n
      go to (108,105,106,107,108),np
  105 nstart = 2
      go to 108
  106 nstart = 2
  107 nstop = np1
  108 nunk = nstop-nstart+1
c
c     define a,b,c coefficients in w-array.
c
      id2 = munk
      id3 = id2+munk
      id4 = id3+munk
      id5 = id4+munk
      id6 = id5+munk
      istart = 1
      a1 = 2./dlrsq
      ij = 0
      if (mbdcnd.eq.3 .or. mbdcnd.eq.4) ij = 1
      if (mbdcnd .le. 4) go to 109
      w(1) = 0.
      w(id2+1) = -2.*a1
      w(id3+1) = 2.*a1
      istart = 2
      ij = 1
  109 do 110 i=istart,munk
         r = a+(i-ij)*deltar
         j = id5+i
         w(j) = r
         j = id6+i
         w(j) = 1./r**2
         w(i) = (r-dlrby2)/(r*dlrsq)
         j = id3+i
         w(j) = (r+dlrby2)/(r*dlrsq)
         k = id6+i
         j = id2+i
         w(j) = -a1+elmbda*w(k)
  110 continue
      go to (114,111,112,113,114,112),mbdcnd
  111 w(id2) = a1
      go to 114
  112 w(id2) = a1
  113 w(id3+1) = a1*istart
  114 continue
c
c     enter boundary data for r-boundaries.
c
      go to (115,115,117,117,119,119),mbdcnd
  115 a1 = w(1)
      do 116 j=nstart,nstop
         f(2,j) = f(2,j)-a1*f(1,j)
  116 continue
      go to 119
  117 a1 = 2.*deltar*w(1)
      do 118 j=nstart,nstop
         f(1,j) = f(1,j)+a1*bda(j)
  118 continue
  119 go to (120,122,122,120,120,122),mbdcnd
  120 a1 = w(id4)
      do 121 j=nstart,nstop
         f(m,j) = f(m,j)-a1*f(mp1,j)
  121 continue
      go to 124
  122 a1 = 2.*deltar*w(id4)
      do 123 j=nstart,nstop
         f(mp1,j) = f(mp1,j)-a1*bdb(j)
  123 continue
c
c     enter boundary data for z-boundaries.
c
  124 a1 = 1./dlthsq
      l = id5-mstart+1
      go to (134,125,125,127,127),np
  125 do 126 i=mstart,mstop
         f(i,2) = f(i,2)-a1*f(i,1)
  126 continue
      go to 129
  127 a1 = 2./deltht
      do 128 i=mstart,mstop
         f(i,1) = f(i,1)+a1*bdc(i)
  128 continue
  129 a1 = 1./dlthsq
      go to (134,130,132,132,130),np
  130 do 131 i=mstart,mstop
         f(i,n) = f(i,n)-a1*f(i,np1)
  131 continue
      go to 134
  132 a1 = 2./deltht
      do 133 i=mstart,mstop
         f(i,np1) = f(i,np1)-a1*bdd(i)
  133 continue
  134 continue
c
c     adjust right side of singular problems to insure existence of a
c     solution.
c
      pertrb = 0.
      if (elmbda) 146,136,135
  135 ierror = 11
      go to 146
  136 w(id5+1) = .5*(w(id5+2)-dlrby2)
      go to (146,146,138,146,146,137),mbdcnd
  137 w(id5+1) = .5*w(id5+1)
  138 go to (140,146,146,139,146),np
  139 a2 = 2.
      go to 141
  140 a2 = 1.
  141 k = id5+munk
      w(k) = .5*(w(k-1)+dlrby2)
      s = 0.
      do 143 i=mstart,mstop
         s1 = 0.
         nsp1 = nstart+1
         nstm1 = nstop-1
         do 142 j=nsp1,nstm1
            s1 = s1+f(i,j)
  142    continue
         k = i+l
         s = s+(a2*s1+f(i,nstart)+f(i,nstop))*w(k)
  143 continue
      s2 = m*a+(.75+(m-1)*(m+1))*dlrby2
      if (mbdcnd .eq. 3) s2 = s2+.25*dlrby2
      s1 = (2.+a2*(nunk-2))*s2
      pertrb = s/s1
      do 145 i=mstart,mstop
         do 144 j=nstart,nstop
            f(i,j) = f(i,j)-pertrb
  144    continue
  145 continue
  146 continue
c
c     multiply i-th equation through by deltht**2 to put equation into
c     correct form for subroutine genbun.
c
      do 148 i=mstart,mstop
         k = i-mstart+1
         w(k) = w(k)*dlthsq
         j = id2+k
         w(j) = w(j)*dlthsq
         j = id3+k
         w(j) = w(j)*dlthsq
         do 147 j=nstart,nstop
            f(i,j) = f(i,j)*dlthsq
  147    continue
  148 continue
      w(1) = 0.
      w(id4) = 0.
c
c     call genbun to solve the system of equations.
c
      call genbun (nbdcnd,nunk,1,munk,w(1),w(id2+1),w(id3+1),idimf,
     1             f(mstart,nstart),ierr1,w(id4+1))
      w(1) = w(id4+1)+3*munk
      if (nbdcnd .ne. 0) go to 150
      do 149 i=mstart,mstop
         f(i,np1) = f(i,1)
  149 continue
  150 continue
      return
      end
