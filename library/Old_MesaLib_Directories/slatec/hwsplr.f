*deck hwsplr
      subroutine hwsplr (a, b, m, mbdcnd, bda, bdb, c, d, n, nbdcnd,
     +   bdc, bdd, elmbda, f, idimf, pertrb, ierror, w)
c***begin prologue  hwsplr
c***purpose  solve a finite difference approximation to the helmholtz
c            equation in polar coordinates.
c***library   slatec (fishpack)
c***category  i2b1a1a
c***type      single precision (hwsplr-s)
c***keywords  elliptic, fishpack, helmholtz, pde, polar
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c     subroutine hwsplr solves a finite difference approximation to the
c     helmholtz equation in polar coordinates:
c
c          (1/r)(d/dr)(r(du/dr)) + (1/r**2)(d/dtheta)(du/dtheta)
c
c                                + lambda*u = f(r,theta).
c
c
c
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
c       indicates the type of boundary condition at r = a and r = b.
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
c            bda(j) = (d/dr)u(a,theta(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bda is a dummy variable.
c
c     bdb
c       a one-dimensional array of length n+1 that specifies the values
c       of the derivative of the solution with respect to r at r = b.
c       when mbdcnd = 2,3, or 6,
c
c            bdb(j) = (d/dr)u(b,theta(j)), j = 1,2,...,n+1  .
c
c       when mbdcnd has any other value, bdb is a dummy variable.
c
c     c,d
c       the range of theta, i.e., c .le. theta .le. d.  c must be less
c       than d.
c
c     n
c       the number of panels into which the interval (c,d) is
c       subdivided.  hence, there will be n+1 grid points in the
c       theta-direction given by theta(j) = c+(j-1)dtheta for
c       j = 1,2,...,n+1, where dtheta = (d-c)/n is the panel width.  n
c       must be greater than 3.
c
c     nbdcnd
c       indicates the type of boundary conditions at theta = c and
c       at theta = d.
c
c       = 0  if the solution is periodic in theta, i.e.,
c            u(i,j) = u(i,n+j).
c       = 1  if the solution is specified at theta = c and theta = d
c            (see note below).
c       = 2  if the solution is specified at theta = c and the
c            derivative of the solution with respect to theta is
c            specified at theta = d (see note below).
c       = 4  if the derivative of the solution with respect to theta is
c            specified at theta = c and the solution is specified at
c            theta = d (see note below).
c
c       note:  when nbdcnd = 1,2, or 4, do not use mbdcnd = 5 or 6
c              (the former indicates that the solution is specified at
c              r = 0, the latter indicates the solution is unspecified
c              at r = 0).  use instead mbdcnd = 1 or 2  .
c
c     bdc
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to theta at
c       theta = c.  when nbdcnd = 3 or 4,
c
c            bdc(i) = (d/dtheta)u(r(i),c), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdc is a dummy variable.
c
c     bdd
c       a one-dimensional array of length m+1 that specifies the values
c       of the derivative of the solution with respect to theta at
c       theta = d.  when nbdcnd = 2 or 3,
c
c            bdd(i) = (d/dtheta)u(r(i),d), i = 1,2,...,m+1  .
c
c       when nbdcnd has any other value, bdd is a dummy variable.
c
c     elmbda
c       the constant lambda in the helmholtz equation.  if
c       lambda .lt. 0, a solution may not exist.  however, hwsplr will
c       attempt to find a solution.
c
c     f
c       a two-dimensional array that specifies the values of the right
c       side of the helmholtz equation and boundary values (if any).
c       for i = 2,3,...,m and j = 2,3,...,n
c
c            f(i,j) = f(r(i),theta(j)).
c
c       on the boundaries f is defined by
c
c            mbdcnd   f(1,j)            f(m+1,j)
c            ------   -------------     -------------
c
c              1      u(a,theta(j))     u(b,theta(j))
c              2      u(a,theta(j))     f(b,theta(j))
c              3      f(a,theta(j))     f(b,theta(j))
c              4      f(a,theta(j))     u(b,theta(j))   j = 1,2,...,n+1
c              5      f(0,0)            u(b,theta(j))
c              6      f(0,0)            f(b,theta(j))
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
c
c     idimf
c       the row (or first) dimension of the array f as it appears in the
c       program calling hwsplr.  this parameter is used to specify the
c       variable dimension of f.  idimf must be at least m+1  .
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space.  w may require up to 4*(n+1) +
c       (13 + int(log2(n+1)))*(m+1) locations.  the actual number of
c       locations used is computed by hwsplr and is returned in location
c       w(1).
c
c
c             * * * * * *   on output     * * * * * *
c
c     f
c       contains the solution u(i,j) of the finite difference
c       approximation for the grid point (r(i),theta(j)),
c       i = 1,2,...,m+1, j = 1,2,...,n+1  .
c
c     pertrb
c       if a combination of periodic, derivative, or unspecified
c       boundary conditions is specified for a poisson equation
c       (lambda = 0), a solution may not exist.  pertrb is a constant,
c       calculated and subtracted from f, which ensures that a solution
c       exists.  hwsplr then computes this solution, which is a least
c       squares solution to the original approximation.  this solution
c       plus any constant is also a solution.  hence, the solution is
c       not unique.  pertrb should be small compared to the right side.
c       otherwise, a solution is obtained to an essentially different
c       problem.  this comparison should always be made to insure that a
c       meaningful solution has been obtained.
c
c     ierror
c       an error flag that indicates invalid input parameters.  except
c       for numbers 0 and 11, a solution is not attempted.
c
c       =  0  no error.
c       =  1  a .lt. 0  .
c       =  2  a .ge. b.
c       =  3  mbdcnd .lt. 1 or mbdcnd .gt. 6  .
c       =  4  c .ge. d.
c       =  5  n .le. 3
c       =  6  nbdcnd .lt. 0 or .gt. 4  .
c       =  7  a = 0, mbdcnd = 3 or 4  .
c       =  8  a .gt. 0, mbdcnd .ge. 5  .
c       =  9  mbdcnd .ge. 5, nbdcnd .ne. 0 and nbdcnd .ne. 3  .
c       = 10  idimf .lt. m+1  .
c       = 11  lambda .gt. 0  .
c       = 12  m .le. 3
c
c       since this is the only means of indicating a possibly incorrect
c       call to hwsplr, the user should test ierror after the call.
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
c     subprograms    hwsplr,genbun,poisd2,poisn2,poisp2,cosgen,merge,
c     required       trix,tri3,pimach
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        standardized april 1, 1973
c                    revised january 1, 1976
c
c     algorithm      the routine defines the finite difference
c                    equations, incorporates boundary data, and adjusts
c                    the right side of singular systems and then calls
c                    genbun to solve the system.
c
c     space          13430(octal) = 5912(decimal)  locations on the ncar
c     required       control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hwsplr is roughly proportional
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
c***end prologue  hwsplr
c
c
      dimension       f(idimf,*)
      dimension       bda(*)     ,bdb(*)     ,bdc(*)     ,bdd(*)     ,
     1                w(*)
c***first executable statement  hwsplr
      ierror = 0
      if (a .lt. 0.) ierror = 1
      if (a .ge. b) ierror = 2
      if (mbdcnd.le.0 .or. mbdcnd.ge.7) ierror = 3
      if (c .ge. d) ierror = 4
      if (n .le. 3) ierror = 5
      if (nbdcnd.le.-1 .or. nbdcnd.ge.5) ierror = 6
      if (a.eq.0. .and. (mbdcnd.eq.3 .or. mbdcnd.eq.4)) ierror = 7
      if (a.gt.0. .and. mbdcnd.ge.5) ierror = 8
      if (mbdcnd.ge.5 .and. nbdcnd.ne.0 .and. nbdcnd.ne.3) ierror = 9
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
      mstop = mp1
      go to (101,105,102,103,104,105),mbdcnd
  101 mstop = m
      go to 105
  102 mstart = 1
      go to 105
  103 mstart = 1
  104 mstop = m
  105 munk = mstop-mstart+1
      nstart = 1
      nstop = n
      go to (109,106,107,108,109),np
  106 nstart = 2
      go to 109
  107 nstart = 2
  108 nstop = np1
  109 nunk = nstop-nstart+1
c
c     define a,b,c coefficients in w-array.
c
      id2 = munk
      id3 = id2+munk
      id4 = id3+munk
      id5 = id4+munk
      id6 = id5+munk
      a1 = 2./dlrsq
      ij = 0
      if (mbdcnd.eq.3 .or. mbdcnd.eq.4) ij = 1
      do 110 i=1,munk
         r = a+(i-ij)*deltar
         j = id5+i
         w(j) = r
         j = id6+i
         w(j) = 1./r**2
         w(i) = (r-dlrby2)/(r*dlrsq)
         j = id3+i
         w(j) = (r+dlrby2)/(r*dlrsq)
         j = id2+i
         w(j) = -a1+elmbda
  110 continue
      go to (114,111,112,113,114,111),mbdcnd
  111 w(id2) = a1
      go to 114
  112 w(id2) = a1
  113 w(id3+1) = a1
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
c     enter boundary data for theta-boundaries.
c
  124 a1 = 1./dlthsq
      l = id5-mstart+1
      lp = id6-mstart+1
      go to (134,125,125,127,127),np
  125 do 126 i=mstart,mstop
         j = i+lp
         f(i,2) = f(i,2)-a1*w(j)*f(i,1)
  126 continue
      go to 129
  127 a1 = 2./deltht
      do 128 i=mstart,mstop
         j = i+lp
         f(i,1) = f(i,1)+a1*w(j)*bdc(i)
  128 continue
  129 a1 = 1./dlthsq
      go to (134,130,132,132,130),np
  130 do 131 i=mstart,mstop
         j = i+lp
         f(i,n) = f(i,n)-a1*w(j)*f(i,np1)
  131 continue
      go to 134
  132 a1 = 2./deltht
      do 133 i=mstart,mstop
         j = i+lp
         f(i,np1) = f(i,np1)-a1*w(j)*bdd(i)
  133 continue
  134 continue
c
c     adjust right side of equation for unknown at pole when have
c     derivative specified boundary conditions.
c
      if (mbdcnd.ge.5 .and. nbdcnd.eq.3)
     1    f(1,1) = f(1,1)-(bdd(2)-bdc(2))*4./(n*deltht*dlrsq)
c
c     adjust right side of singular problems to insure existence of a
c     solution.
c
      pertrb = 0.
      if (elmbda) 144,136,135
  135 ierror = 11
      go to 144
  136 if (nbdcnd.ne.0 .and. nbdcnd.ne.3) go to 144
      s2 = 0.
      go to (144,144,137,144,144,138),mbdcnd
  137 w(id5+1) = .5*(w(id5+2)-dlrby2)
      s2 = .25*deltar
  138 a2 = 2.
      if (nbdcnd .eq. 0) a2 = 1.
      j = id5+munk
      w(j) = .5*(w(j-1)+dlrby2)
      s = 0.
      do 140 i=mstart,mstop
         s1 = 0.
         ij = nstart+1
         k = nstop-1
         do 139 j=ij,k
            s1 = s1+f(i,j)
  139    continue
         j = i+l
         s = s+(a2*s1+f(i,nstart)+f(i,nstop))*w(j)
  140 continue
      s2 = m*a+deltar*((m-1)*(m+1)*.5+.25)+s2
      s1 = (2.+a2*(nunk-2))*s2
      if (mbdcnd .eq. 3) go to 141
      s2 = n*a2*deltar/8.
      s = s+f(1,1)*s2
      s1 = s1+s2
  141 continue
      pertrb = s/s1
      do 143 i=mstart,mstop
         do 142 j=nstart,nstop
            f(i,j) = f(i,j)-pertrb
  142    continue
  143 continue
  144 continue
c
c     multiply i-th equation through by (r(i)*deltht)**2.
c
      do 146 i=mstart,mstop
         k = i-mstart+1
         j = i+lp
         a1 = dlthsq/w(j)
         w(k) = a1*w(k)
         j = id2+k
         w(j) = a1*w(j)
         j = id3+k
         w(j) = a1*w(j)
         do 145 j=nstart,nstop
            f(i,j) = a1*f(i,j)
  145    continue
  146 continue
      w(1) = 0.
      w(id4) = 0.
c
c     call genbun to solve the system of equations.
c
      call genbun (nbdcnd,nunk,1,munk,w(1),w(id2+1),w(id3+1),idimf,
     1             f(mstart,nstart),ierr1,w(id4+1))
      iwstor = w(id4+1)+3*munk
      go to (157,157,157,157,148,147),mbdcnd
c
c     adjust the solution as necessary for the problems where a = 0.
c
  147 if (elmbda .ne. 0.) go to 148
      ypole = 0.
      go to 155
  148 continue
      j = id5+munk
      w(j) = w(id2)/w(id3)
      do 149 ip=3,munk
         i = munk-ip+2
         j = id5+i
         lp = id2+i
         k = id3+i
         w(j) = w(i)/(w(lp)-w(k)*w(j+1))
  149 continue
      w(id5+1) = -.5*dlthsq/(w(id2+1)-w(id3+1)*w(id5+2))
      do 150 i=2,munk
         j = id5+i
         w(j) = -w(j)*w(j-1)
  150 continue
      s = 0.
      do 151 j=nstart,nstop
         s = s+f(2,j)
  151 continue
      a2 = nunk
      if (nbdcnd .eq. 0) go to 152
      s = s-.5*(f(2,nstart)+f(2,nstop))
      a2 = a2-1.
  152 ypole = (.25*dlrsq*f(1,1)-s/a2)/(w(id5+1)-1.+elmbda*dlrsq*.25)
      do 154 i=mstart,mstop
         k = l+i
         do 153 j=nstart,nstop
            f(i,j) = f(i,j)+ypole*w(k)
  153    continue
  154 continue
  155 do 156 j=1,np1
         f(1,j) = ypole
  156 continue
  157 continue
      if (nbdcnd .ne. 0) go to 159
      do 158 i=mstart,mstop
         f(i,np1) = f(i,1)
  158 continue
  159 continue
      w(1) = iwstor
      return
      end
