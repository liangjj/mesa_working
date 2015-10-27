*deck genbun
      subroutine genbun (nperod, n, mperod, m, a, b, c, idimy, y,
     +   ierror, w)
c***begin prologue  genbun
c***purpose  solve by a cyclic reduction algorithm the linear system
c            of equations that results from a finite difference
c            approximation to certain 2-d elliptic pde's on a centered
c            grid .
c***library   slatec (fishpack)
c***category  i2b4b
c***type      single precision (genbun-s, cmgnbn-c)
c***keywords  elliptic, fishpack, pde, tridiagonal
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c     subroutine genbun solves the linear system of equations
c
c          a(i)*x(i-1,j) + b(i)*x(i,j) + c(i)*x(i+1,j)
c
c          + x(i,j-1) - 2.*x(i,j) + x(i,j+1) = y(i,j)
c
c               for i = 1,2,...,m  and  j = 1,2,...,n.
c
c     the indices i+1 and i-1 are evaluated modulo m, i.e.,
c     x(0,j) = x(m,j) and x(m+1,j) = x(1,j), and x(i,0) may be equal to
c     0, x(i,2), or x(i,n) and x(i,n+1) may be equal to 0, x(i,n-1), or
c     x(i,1) depending on an input parameter.
c
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c     nperod
c       indicates the values that x(i,0) and x(i,n+1) are assumed to
c       have.
c
c       = 0  if x(i,0) = x(i,n) and x(i,n+1) = x(i,1).
c       = 1  if x(i,0) = x(i,n+1) = 0  .
c       = 2  if x(i,0) = 0 and x(i,n+1) = x(i,n-1).
c       = 3  if x(i,0) = x(i,2) and x(i,n+1) = x(i,n-1).
c       = 4  if x(i,0) = x(i,2) and x(i,n+1) = 0.
c
c     n
c       the number of unknowns in the j-direction.  n must be greater
c       than 2.
c
c     mperod
c       = 0 if a(1) and c(m) are not zero.
c       = 1 if a(1) = c(m) = 0.
c
c     m
c       the number of unknowns in the i-direction.  m must be greater
c       than 2.
c
c     a,b,c
c       one-dimensional arrays of length m that specify the
c       coefficients in the linear equations given above.  if mperod = 0
c       the array elements must not depend upon the index i, but must be
c       constant.  specifically, the subroutine checks the following
c       condition
c
c             a(i) = c(1)
c             c(i) = c(1)
c             b(i) = b(1)
c
c       for i=1,2,...,m.
c
c     idimy
c       the row (or first) dimension of the two-dimensional array y as
c       it appears in the program calling genbun.  this parameter is
c       used to specify the variable dimension of y.  idimy must be at
c       least m.
c
c     y
c       a two-dimensional array that specifies the values of the right
c       side of the linear system of equations given above.  y must be
c       dimensioned at least m*n.
c
c     w
c       a one-dimensional array that must be provided by the user for
c       work space.  w may require up to 4*n + (10 + int(log2(n)))*m
c       locations.  the actual number of locations used is computed by
c       genbun and is returned in location w(1).
c
c
c             * * * * * *   on output     * * * * * *
c
c     y
c       contains the solution x.
c
c     ierror
c       an error flag that indicates invalid input parameters.  except
c       for number zero, a solution is not attempted.
c
c       = 0  no error.
c       = 1  m .le. 2
c       = 2  n .le. 2
c       = 3  idimy .lt. m
c       = 4  nperod .lt. 0 or nperod .gt. 4
c       = 5  mperod .lt. 0 or mperod .gt. 1
c       = 6  a(i) .ne. c(1) or c(i) .ne. c(1) or b(i) .ne. b(1) for
c            some i=1,2,...,m.
c       = 7  a(1) .ne. 0 or c(m) .ne. 0 and mperod = 1
c
c     w
c       w(1) contains the required length of w.
c
c *long description:
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   a(m),b(m),c(m),y(idimy,n),w(see parameter list)
c     arguments
c
c     latest         june 1, 1976
c     revision
c
c     subprograms    genbun,poisd2,poisn2,poisp2,cosgen,merge,trix,tri3,
c     required       pimach
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
c     history        standardized april 1, 1973
c                    revised august 20,1973
c                    revised january 1, 1976
c
c     algorithm      the linear system is solved by a cyclic reduction
c                    algorithm described in the reference.
c
c     space          4944(decimal) = 11520(octal) locations on the ncar
c     required       control data 7600.
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine genbun is roughly proportional
c                    to m*n*log2(n), but also depends on the input
c                    parameter nperod.  some typical values are listed
c                    in the table below.  more comprehensive timing
c                    charts may be found in the reference.
c                       to measure the accuracy of the algorithm a
c                    uniform random number generator was used to create
c                    a solution array x for the system given in the
c                    'purpose' with
c
c                       a(i) = c(i) = -0.5*b(i) = 1,       i=1,2,...,m
c
c                    and, when mperod = 1
c
c                       a(1) = c(m) = 0
c                       a(m) = c(1) = 2.
c
c                    the solution x was substituted into the given sys-
c                    tem and, using double precision, a right side y was
c                    computed.  using this array y subroutine genbun was
c                    called to produce an approximate solution z.  then
c                    the relative error, defined as
c
c                       e = max(abs(z(i,j)-x(i,j)))/max(abs(x(i,j)))
c
c                    where the two maxima are taken over all i=1,2,...,m
c                    and j=1,2,...,n, was computed.  the value of e is
c                    given in the table below for some typical values of
c                    m and n.
c
c
c                       m (=n)    mperod    nperod    t(msecs)    e
c                       ------    ------    ------    --------  ------
c
c                         31        0         0          36     6.e-14
c                         31        1         1          21     4.e-13
c                         31        1         3          41     3.e-13
c                         32        0         0          29     9.e-14
c                         32        1         1          32     3.e-13
c                         32        1         3          48     1.e-13
c                         33        0         0          36     9.e-14
c                         33        1         1          30     4.e-13
c                         33        1         3          34     1.e-13
c                         63        0         0         150     1.e-13
c                         63        1         1          91     1.e-12
c                         63        1         3         173     2.e-13
c                         64        0         0         122     1.e-13
c                         64        1         1         128     1.e-12
c                         64        1         3         199     6.e-13
c                         65        0         0         143     2.e-13
c                         65        1         1         120     1.e-12
c                         65        1         3         138     4.e-13
c
c     portability    american national standards institute fortran.
c                    the machine dependent constant pi is defined in
c                    function pimach.
c
c     required       cos
c     resident
c     routines
c
c     reference      sweet, r., 'a cyclic reduction algorithm for
c                    solving block tridiagonal systems of arbitrary
c                    dimensions,' siam j. on numer. anal.,
c                    14(sept., 1977), pp. 706-720.
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c***references  r. sweet, a cyclic reduction algorithm for solving
c                 block tridiagonal systems of arbitrary dimensions,
c                 siam journal on numerical analysis 14, (september
c                 1977), pp. 706-720.
c***routines called  poisd2, poisn2, poisp2
c***revision history  (yymmdd)
c   801001  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  genbun
c
c
      dimension       y(idimy,*)
      dimension       w(*)       ,b(*)       ,a(*)       ,c(*)
c***first executable statement  genbun
      ierror = 0
      if (m .le. 2) ierror = 1
      if (n .le. 2) ierror = 2
      if (idimy .lt. m) ierror = 3
      if (nperod.lt.0 .or. nperod.gt.4) ierror = 4
      if (mperod.lt.0 .or. mperod.gt.1) ierror = 5
      if (mperod .eq. 1) go to 102
      do 101 i=2,m
         if (a(i) .ne. c(1)) go to 103
         if (c(i) .ne. c(1)) go to 103
         if (b(i) .ne. b(1)) go to 103
  101 continue
      go to 104
  102 if (a(1).ne.0. .or. c(m).ne.0.) ierror = 7
      go to 104
  103 ierror = 6
  104 if (ierror .ne. 0) return
      mp1 = m+1
      iwba = mp1
      iwbb = iwba+m
      iwbc = iwbb+m
      iwb2 = iwbc+m
      iwb3 = iwb2+m
      iww1 = iwb3+m
      iww2 = iww1+m
      iww3 = iww2+m
      iwd = iww3+m
      iwtcos = iwd+m
      iwp = iwtcos+4*n
      do 106 i=1,m
         k = iwba+i-1
         w(k) = -a(i)
         k = iwbc+i-1
         w(k) = -c(i)
         k = iwbb+i-1
         w(k) = 2.-b(i)
         do 105 j=1,n
            y(i,j) = -y(i,j)
  105    continue
  106 continue
      mp = mperod+1
      np = nperod+1
      go to (114,107),mp
  107 go to (108,109,110,111,123),np
  108 call poisp2 (m,n,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2),
     1             w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos),
     2             w(iwp))
      go to 112
  109 call poisd2 (m,n,1,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iww1),
     1             w(iwd),w(iwtcos),w(iwp))
      go to 112
  110 call poisn2 (m,n,1,2,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2),
     1             w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos),
     2             w(iwp))
      go to 112
  111 call poisn2 (m,n,1,1,w(iwba),w(iwbb),w(iwbc),y,idimy,w,w(iwb2),
     1             w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos),
     2             w(iwp))
  112 ipstor = w(iww1)
      irev = 2
      if (nperod .eq. 4) go to 124
  113 go to (127,133),mp
  114 continue
c
c     reorder unknowns when mp =0
c
      mh = (m+1)/2
      mhm1 = mh-1
      modd = 1
      if (mh*2 .eq. m) modd = 2
      do 119 j=1,n
         do 115 i=1,mhm1
            mhpi = mh+i
            mhmi = mh-i
            w(i) = y(mhmi,j)-y(mhpi,j)
            w(mhpi) = y(mhmi,j)+y(mhpi,j)
  115    continue
         w(mh) = 2.*y(mh,j)
         go to (117,116),modd
  116    w(m) = 2.*y(m,j)
  117    continue
         do 118 i=1,m
            y(i,j) = w(i)
  118    continue
  119 continue
      k = iwbc+mhm1-1
      i = iwba+mhm1
      w(k) = 0.
      w(i) = 0.
      w(k+1) = 2.*w(k+1)
      go to (120,121),modd
  120 continue
      k = iwbb+mhm1-1
      w(k) = w(k)-w(i-1)
      w(iwbc-1) = w(iwbc-1)+w(iwbb-1)
      go to 122
  121 w(iwbb-1) = w(k+1)
  122 continue
      go to 107
c
c     reverse columns when nperod = 4.
c
  123 irev = 1
      nby2 = n/2
  124 do 126 j=1,nby2
         mskip = n+1-j
         do 125 i=1,m
            a1 = y(i,j)
            y(i,j) = y(i,mskip)
            y(i,mskip) = a1
  125    continue
  126 continue
      go to (110,113),irev
  127 continue
      do 132 j=1,n
         do 128 i=1,mhm1
            mhmi = mh-i
            mhpi = mh+i
            w(mhmi) = .5*(y(mhpi,j)+y(i,j))
            w(mhpi) = .5*(y(mhpi,j)-y(i,j))
  128    continue
         w(mh) = .5*y(mh,j)
         go to (130,129),modd
  129    w(m) = .5*y(m,j)
  130    continue
         do 131 i=1,m
            y(i,j) = w(i)
  131    continue
  132 continue
  133 continue
c
c     return storage requirements for w array.
c
      w(1) = ipstor+iwp-1
      return
      end
