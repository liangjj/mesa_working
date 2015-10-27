*deck poistg
      subroutine poistg (nperod, n, mperod, m, a, b, c, idimy, y,
     +   ierror, w)
c***begin prologue  poistg
c***purpose  solve a block tridiagonal system of linear equations
c            that results from a staggered grid finite difference
c            approximation to 2-d elliptic pde's.
c***library   slatec (fishpack)
c***category  i2b4b
c***type      single precision (poistg-s)
c***keywords  elliptic, fishpack, helmholtz, pde, tridiagonal
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c     subroutine poistg solves the linear system of equations
c
c       a(i)*x(i-1,j) + b(i)*x(i,j) + c(i)*x(i+1,j)
c       + x(i,j-1) - 2.*x(i,j) + x(i,j+1) = y(i,j)
c
c       for i=1,2,...,m and j=1,2,...,n.
c
c     the indices i+1 and i-1 are evaluated modulo m, i.e.
c     x(0,j) = x(m,j) and x(m+1,j) = x(1,j), and x(i,0) may be equal to
c     x(i,1) or -x(i,1) and x(i,n+1) may be equal to x(i,n) or -x(i,n)
c     depending on an input parameter.
c
c
c     * * * * * * * *    parameter description     * * * * * * * * * *
c
c             * * * * * *   on input    * * * * * *
c
c   nperod
c     indicates the values which x(i,0) and x(i,n+1) are assumed
c     to have.
c     = 1 if x(i,0) = -x(i,1) and x(i,n+1) = -x(i,n)
c     = 2 if x(i,0) = -x(i,1) and x(i,n+1) =  x(i,n)
c     = 3 if x(i,0) =  x(i,1) and x(i,n+1) =  x(i,n)
c     = 4 if x(i,0) =  x(i,1) and x(i,n+1) = -x(i,n)
c
c   n
c     the number of unknowns in the j-direction.  n must
c     be greater than 2.
c
c   mperod
c     = 0 if a(1) and c(m) are not zero
c     = 1 if a(1) = c(m) = 0
c
c   m
c     the number of unknowns in the i-direction.  m must
c     be greater than 2.
c
c   a,b,c
c     one-dimensional arrays of length m that specify the coefficients
c     in the linear equations given above.  if mperod = 0 the array
c     elements must not depend on the index i, but must be constant.
c     specifically, the subroutine checks the following condition
c
c           a(i) = c(1)
c           b(i) = b(1)
c           c(i) = c(1)
c
c     for i = 1, 2, ..., m.
c
c   idimy
c     the row (or first) dimension of the two-dimensional array y as
c     it appears in the program calling poistg.  this parameter is
c     used to specify the variable dimension of y.  idimy must be at
c     least m.
c
c   y
c     a two-dimensional array that specifies the values of the
c     right side of the linear system of equations given above.
c     y must be dimensioned at least m x n.
c
c   w
c     a one-dimensional work array that must be provided by the user
c     for work space.  w may require up to 9m + 4n + m(int(log2(n)))
c     locations.  the actual number of locations used is computed by
c     poistg and returned in location w(1).
c
c
c             * * * * * *   on output     * * * * * *
c
c   y
c     contains the solution x.
c
c   ierror
c     an error flag that indicates invalid input parameters.  except
c     for number zero, a solution is not attempted.
c     = 0  no error
c     = 1  if m .le. 2
c     = 2  if n .le. 2
c     = 3  idimy .lt. m
c     = 4  if nperod .lt. 1 or nperod .gt. 4
c     = 5  if mperod .lt. 0 or mperod .gt. 1
c     = 6  if mperod = 0 and
c          a(i) .ne. c(1) or b(i) .ne. b(1) or c(i) .ne. c(1)
c          for some i = 1, 2, ..., m.
c       = 7 if mperod .eq. 1 .and. (a(1).ne.0 .or. c(m).ne.0)
c
c   w
c     w(1) contains the required length of w.
c
c *long description:
c
c     * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   a(m),b(m),c(m),y(idimy,n),
c     arguments      w(see argument list)
c
c     latest         june 1, 1977
c     revision
c
c     subprograms    poistg,postg2,cosgen,merge,trix,tri3,pimach
c     required
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
c     history        written by roland sweet in 1973
c                    revised by roland sweet in 1977
c
c
c     space          3297(decimal) = 6341(octal) locations on the
c     required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine poistg is roughly proportional
c                    to m*n*log2(n).  some typical values are listed
c                    in the table below.  more comprehensive timing
c                    charts may be found in the reference.
c                       to measure the accuracy of the algorithm a
c                    uniform random number generator was used to create
c                    a solution array x for the system given in the
c                    'purpose ' with
c
c                       a(i) = c(i) = -0.5*b(i) = 1,       i=1,2,...,m
c
c                    and, when mperod = 1
c
c                       a(1) = c(m) = 0
c                       b(1) = b(m) =-1.
c
c                    the solution x was substituted into the given sys-
c                    tem and, using double precision, a right side y was
c                    computed.  using this array y subroutine poistg was
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
c                         31        0-1       1-4        45     9.e-13
c                         31        1         1          21     4.e-13
c                         31        1         3          41     3.e-13
c                         32        0-1       1-4        51     3.e-12
c                         32        1         1          32     3.e-13
c                         32        1         3          48     1.e-13
c                         33        0-1       1-4        42     1.e-12
c                         33        1         1          30     4.e-13
c                         33        1         3          34     1.e-13
c                         63        0-1       1-4       186     3.e-12
c                         63        1         1          91     1.e-12
c                         63        1         3         173     2.e-13
c                         64        0-1       1-4       209     4.e-12
c                         64        1         1         128     1.e-12
c                         64        1         3         199     6.e-13
c                         65        0-1       1-4       143     2.e-13
c                         65        1         1         160     1.e-11
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
c***routines called  postg2
c***revision history  (yymmdd)
c   801001  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  poistg
c
c
      dimension       y(idimy,*)
      dimension       w(*)       ,b(*)       ,a(*)       ,c(*)
c***first executable statement  poistg
      ierror = 0
      if (m .le. 2) ierror = 1
      if (n .le. 2) ierror = 2
      if (idimy .lt. m) ierror = 3
      if (nperod.lt.1 .or. nperod.gt.4) ierror = 4
      if (mperod.lt.0 .or. mperod.gt.1) ierror = 5
      if (mperod .eq. 1) go to 103
      do 101 i=1,m
         if (a(i) .ne. c(1)) go to 102
         if (c(i) .ne. c(1)) go to 102
         if (b(i) .ne. b(1)) go to 102
  101 continue
      go to 104
  102 ierror = 6
      return
  103 if (a(1).ne.0. .or. c(m).ne.0.) ierror = 7
  104 if (ierror .ne. 0) return
      iwba = m+1
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
      np = nperod
      mp = mperod+1
      go to (110,107),mp
  107 continue
      go to (108,108,108,119),nperod
  108 continue
      call postg2 (np,n,m,w(iwba),w(iwbb),w(iwbc),idimy,y,w,w(iwb2),
     1             w(iwb3),w(iww1),w(iww2),w(iww3),w(iwd),w(iwtcos),
     2             w(iwp))
      ipstor = w(iww1)
      irev = 2
      if (nperod .eq. 4) go to 120
  109 continue
      go to (123,129),mp
  110 continue
c
c     reorder unknowns when mp =0
c
      mh = (m+1)/2
      mhm1 = mh-1
      modd = 1
      if (mh*2 .eq. m) modd = 2
      do 115 j=1,n
         do 111 i=1,mhm1
            mhpi = mh+i
            mhmi = mh-i
            w(i) = y(mhmi,j)-y(mhpi,j)
            w(mhpi) = y(mhmi,j)+y(mhpi,j)
  111    continue
         w(mh) = 2.*y(mh,j)
         go to (113,112),modd
  112    w(m) = 2.*y(m,j)
  113    continue
         do 114 i=1,m
            y(i,j) = w(i)
  114    continue
  115 continue
      k = iwbc+mhm1-1
      i = iwba+mhm1
      w(k) = 0.
      w(i) = 0.
      w(k+1) = 2.*w(k+1)
      go to (116,117),modd
  116 continue
      k = iwbb+mhm1-1
      w(k) = w(k)-w(i-1)
      w(iwbc-1) = w(iwbc-1)+w(iwbb-1)
      go to 118
  117 w(iwbb-1) = w(k+1)
  118 continue
      go to 107
  119 continue
c
c     reverse columns when nperod = 4.
c
      irev = 1
      nby2 = n/2
      np = 2
  120 do 122 j=1,nby2
         mskip = n+1-j
         do 121 i=1,m
            a1 = y(i,j)
            y(i,j) = y(i,mskip)
            y(i,mskip) = a1
  121    continue
  122 continue
      go to (108,109),irev
  123 continue
      do 128 j=1,n
         do 124 i=1,mhm1
            mhmi = mh-i
            mhpi = mh+i
            w(mhmi) = .5*(y(mhpi,j)+y(i,j))
            w(mhpi) = .5*(y(mhpi,j)-y(i,j))
  124    continue
         w(mh) = .5*y(mh,j)
         go to (126,125),modd
  125    w(m) = .5*y(m,j)
  126    continue
         do 127 i=1,m
            y(i,j) = w(i)
  127    continue
  128 continue
  129 continue
c
c     return storage requirements for w array.
c
      w(1) = ipstor+iwp-1
      return
      end
