*deck pois3d
      subroutine pois3d (lperod, l, c1, mperod, m, c2, nperod, n, a, b,
     +   c, ldimf, mdimf, f, ierror, w)
c***begin prologue  pois3d
c***purpose  solve a three-dimensional block tridiagonal linear system
c            which arises from a finite difference approximation to a
c            three-dimensional poisson equation using the fourier
c            transform package fftpak written by paul swarztrauber.
c***library   slatec (fishpack)
c***category  i2b4b
c***type      single precision (pois3d-s)
c***keywords  elliptic pde, fishpack, helmholtz, poisson
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c     subroutine pois3d solves the linear system of equations
c
c       c1*(x(i-1,j,k)-2.*x(i,j,k)+x(i+1,j,k))
c     + c2*(x(i,j-1,k)-2.*x(i,j,k)+x(i,j+1,k))
c     + a(k)*x(i,j,k-1)+b(k)*x(i,j,k)+c(k)*x(i,j,k+1) = f(i,j,k)
c
c     for  i=1,2,...,l , j=1,2,...,m , and k=1,2,...,n .
c
c     the indices k-1 and k+1 are evaluated modulo n, i.e.
c     x(i,j,0) = x(i,j,n) and x(i,j,n+1) = x(i,j,1). the unknowns
c     x(0,j,k), x(l+1,j,k), x(i,0,k), and x(i,m+1,k) are assumed to take
c     on certain prescribed values described below.
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c    * * * * * * * *    parameter description     * * * * * * * * * *
c
c
c            * * * * * *   on input    * * * * * *
c
c     lperod   indicates the values that x(0,j,k) and x(l+1,j,k) are
c              assumed to have.
c
c              = 0  if x(0,j,k) = x(l,j,k) and x(l+1,j,k) = x(1,j,k).
c              = 1  if x(0,j,k) = x(l+1,j,k) = 0.
c              = 2  if x(0,j,k) = 0  and x(l+1,j,k) = x(l-1,j,k).
c              = 3  if x(0,j,k) = x(2,j,k) and x(l+1,j,k) = x(l-1,j,k).
c              = 4  if x(0,j,k) = x(2,j,k) and x(l+1,j,k) = 0.
c
c     l        the number of unknowns in the i-direction. l must be at
c              least 3.
c
c     c1       the real constant that appears in the above equation.
c
c     mperod   indicates the values that x(i,0,k) and x(i,m+1,k) are
c              assumed to have.
c
c              = 0  if x(i,0,k) = x(i,m,k) and x(i,m+1,k) = x(i,1,k).
c              = 1  if x(i,0,k) = x(i,m+1,k) = 0.
c              = 2  if x(i,0,k) = 0 and x(i,m+1,k) = x(i,m-1,k).
c              = 3  if x(i,0,k) = x(i,2,k) and x(i,m+1,k) = x(i,m-1,k).
c              = 4  if x(i,0,k) = x(i,2,k) and x(i,m+1,k) = 0.
c
c     m        the number of unknowns in the j-direction. m must be at
c              least 3.
c
c     c2       the real constant which appears in the above equation.
c
c     nperod   = 0  if a(1) and c(n) are not zero.
c              = 1  if a(1) = c(n) = 0.
c
c     n        the number of unknowns in the k-direction. n must be at
c              least 3.
c
c
c     a,b,c    one-dimensional arrays of length n that specify the
c              coefficients in the linear equations given above.
c
c              if nperod = 0 the array elements must not depend upon the
c              index k, but must be constant.  specifically, the
c              subroutine checks the following condition
c
c                          a(k) = c(1)
c                          c(k) = c(1)
c                          b(k) = b(1)
c
c                  for k=1,2,...,n.
c
c     ldimf    the row (or first) dimension of the three-dimensional
c              array f as it appears in the program calling pois3d.
c              this parameter is used to specify the variable dimension
c              of f.  ldimf must be at least l.
c
c     mdimf    the column (or second) dimension of the three-dimensional
c              array f as it appears in the program calling pois3d.
c              this parameter is used to specify the variable dimension
c              of f.  mdimf must be at least m.
c
c     f        a three-dimensional array that specifies the values of
c              the right side of the linear system of equations given
c              above.  f must be dimensioned at least l x m x n.
c
c     w        a one-dimensional array that must be provided by the
c              user for work space.  the length of w must be at least
c              30 + l + m + 2*n + max(l,m,n) +
c              7*(int((l+1)/2) + int((m+1)/2)).
c
c
c            * * * * * *   on output   * * * * * *
c
c     f        contains the solution x.
c
c     ierror   an error flag that indicates invalid input parameters.
c              except for number zero, a solution is not attempted.
c              = 0  no error
c              = 1  if lperod .lt. 0 or .gt. 4
c              = 2  if l .lt. 3
c              = 3  if mperod .lt. 0 or .gt. 4
c              = 4  if m .lt. 3
c              = 5  if nperod .lt. 0 or .gt. 1
c              = 6  if n .lt. 3
c              = 7  if ldimf .lt. l
c              = 8  if mdimf .lt. m
c              = 9  if a(k) .ne. c(1) or c(k) .ne. c(1) or b(i) .ne.b(1)
c                      for some k=1,2,...,n.
c              = 10 if nperod = 1 and a(1) .ne. 0 or c(n) .ne. 0
c
c              since this is the only means of indicating a possibly
c              incorrect call to pois3d, the user should test ierror
c              after the call.
c
c *long description:
c
c    * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   a(n),b(n),c(n),f(ldimf,mdimf,n),
c     arguments      w(see argument list)
c
c     latest         december 1, 1978
c     revision
c
c     subprograms    pois3d,pos3d1,tridq,rffti,rfftf,rfftf1,rfftb,
c     required       rfftb1,costi,cost,sinti,sint,cosqi,cosqf,cosqf1
c                    cosqb,cosqb1,sinqi,sinqf,sinqb,cffti,cffti1,
c                    cfftb,cfftb1,passb2,passb3,passb4,passb,cfftf,
c                    cfftf1,passf1,passf2,passf3,passf4,passf,pimach,
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
c     history        written by roland sweet at ncar in july 1977
c
c     algorithm      this subroutine solves three-dimensional block
c                    tridiagonal linear systems arising from finite
c                    difference approximations to three-dimensional
c                    poisson equations using the fourier transform
c                    package fftpak written by paul swarztrauber.
c
c     space          6561(decimal) = 14641(octal) locations on the
c     required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine pois3d is roughly proportional
c                    to l*m*n*(log2(l)+log2(m)+5), but also depends on
c                    input parameters lperod and mperod.  some typical
c                    values are listed in the table below when nperod=0.
c                       to measure the accuracy of the algorithm a
c                    uniform random number generator was used to create
c                    a solution array x for the system given in the
c                    'purpose' with
c
c                       a(k) = c(k) = -0.5*b(k) = 1,       k=1,2,...,n
c
c                    and, when nperod = 1
c
c                       a(1) = c(n) = 0
c                       a(n) = c(1) = 2.
c
c                    the solution x was substituted into the given sys-
c                    tem and, using double precision, a right side y was
c                    computed.  using this array y subroutine pois3d was
c                    called to produce an approximate solution z.  then
c                    the relative error, defined as
c
c                    e = max(abs(z(i,j,k)-x(i,j,k)))/max(abs(x(i,j,k)))
c
c                    where the two maxima are taken over i=1,2,...,l,
c                    j=1,2,...,m and k=1,2,...,n, was computed.  the
c                    value of e is given in the table below for some
c                    typical values of l,m and n.
c
c
c                       l(=m=n)   lperod    mperod    t(msecs)    e
c                       ------    ------    ------    --------  ------
c
c                         16        0         0         272     1.e-13
c                         15        1         1         287     4.e-13
c                         17        3         3         338     2.e-13
c                         32        0         0        1755     2.e-13
c                         31        1         1        1894     2.e-12
c                         33        3         3        2042     7.e-13
c
c
c     portability    american national standards institute fortran.
c                    the machine dependent constant pi is defined in
c                    function pimach.
c
c     required       cos,sin,atan
c     resident
c     routines
c
c     reference      none
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c***references  (none)
c***routines called  pos3d1
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  pois3d
      dimension       a(*)       ,b(*)       ,c(*)       ,
     1                f(ldimf,mdimf,*)       ,w(*)       ,save(6)
c***first executable statement  pois3d
      lp = lperod+1
      mp = mperod+1
      np = nperod+1
c
c     check for invalid input.
c
      ierror = 0
      if (lp.lt.1 .or. lp.gt.5) ierror = 1
      if (l .lt. 3) ierror = 2
      if (mp.lt.1 .or. mp.gt.5) ierror = 3
      if (m .lt. 3) ierror = 4
      if (np.lt.1 .or. np.gt.2) ierror = 5
      if (n .lt. 3) ierror = 6
      if (ldimf .lt. l) ierror = 7
      if (mdimf .lt. m) ierror = 8
      if (np .ne. 1) go to 103
      do 101 k=1,n
         if (a(k) .ne. c(1)) go to 102
         if (c(k) .ne. c(1)) go to 102
         if (b(k) .ne. b(1)) go to 102
  101 continue
      go to 104
  102 ierror = 9
  103 if (nperod.eq.1 .and. (a(1).ne.0. .or. c(n).ne.0.)) ierror = 10
  104 if (ierror .ne. 0) go to 122
      iwyrt = l+1
      iwt = iwyrt+m
      iwd = iwt+max(l,m,n)+1
      iwbb = iwd+n
      iwx = iwbb+n
      iwy = iwx+7*((l+1)/2)+15
      go to (105,114),np
c
c     reorder unknowns when nperod = 0.
c
  105 nh = (n+1)/2
      nhm1 = nh-1
      nodd = 1
      if (2*nh .eq. n) nodd = 2
      do 111 i=1,l
         do 110 j=1,m
            do 106 k=1,nhm1
               nhpk = nh+k
               nhmk = nh-k
               w(k) = f(i,j,nhmk)-f(i,j,nhpk)
               w(nhpk) = f(i,j,nhmk)+f(i,j,nhpk)
  106       continue
            w(nh) = 2.*f(i,j,nh)
            go to (108,107),nodd
  107       w(n) = 2.*f(i,j,n)
  108       do 109 k=1,n
               f(i,j,k) = w(k)
  109       continue
  110    continue
  111 continue
      save(1) = c(nhm1)
      save(2) = a(nh)
      save(3) = c(nh)
      save(4) = b(nhm1)
      save(5) = b(n)
      save(6) = a(n)
      c(nhm1) = 0.
      a(nh) = 0.
      c(nh) = 2.*c(nh)
      go to (112,113),nodd
  112 b(nhm1) = b(nhm1)-a(nh-1)
      b(n) = b(n)+a(n)
      go to 114
  113 a(n) = c(nh)
  114 continue
      call pos3d1 (lp,l,mp,m,n,a,b,c,ldimf,mdimf,f,w,w(iwyrt),w(iwt),
     1             w(iwd),w(iwx),w(iwy),c1,c2,w(iwbb))
      go to (115,122),np
  115 do 121 i=1,l
         do 120 j=1,m
            do 116 k=1,nhm1
               nhmk = nh-k
               nhpk = nh+k
               w(nhmk) = .5*(f(i,j,nhpk)+f(i,j,k))
               w(nhpk) = .5*(f(i,j,nhpk)-f(i,j,k))
  116       continue
            w(nh) = .5*f(i,j,nh)
            go to (118,117),nodd
  117       w(n) = .5*f(i,j,n)
  118       do 119 k=1,n
               f(i,j,k) = w(k)
  119       continue
  120    continue
  121 continue
      c(nhm1) = save(1)
      a(nh) = save(2)
      c(nh) = save(3)
      b(nhm1) = save(4)
      b(n) = save(5)
      a(n) = save(6)
  122 continue
      return
      end
