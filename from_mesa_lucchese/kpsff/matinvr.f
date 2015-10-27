      subroutine matinvr ( a, n, b, m, det, nmax )
      implicit real*8(a-h,o-z)
      real*8      a, b, det, swap, t
      real*8  amax, temp
      real*8  pivot
      dimension a(nmax,n), b(nmax,m)
      common / f402 / pivot(4096), index(4096)
c     -----------------------------------------------------------------
c     initialize determinant and pivot element array
c     -----------------------------------------------------------------
      det = 1e0
      do 20 i = 1, n
      pivot(i) = 0.
   20 continue
c     perform successive pivot operations ( grand loop )
c     -----------------------------------------------------------------
      do 550 i = 1, n
c     -----------------------------------------------------------------
c     search for pivot element and extend determinant partial product
c     -----------------------------------------------------------------
      amax = 0e0
      do 105 j = 1, n
      if ( pivot(j) .ne. 0. ) go to 105
      do 100 k = 1, n
      if ( pivot(k) .ne. 0. ) go to 100
      temp = abs ( a(j,k) )
      if ( temp .lt. amax ) go to 100
      ir = j
      ic = k
      amax = temp
  100 continue
  105 continue
      index(i) = 4096*ir + ic
      j = ir
      t = a(j,ic)
      det = det*t
c     -----------------------------------------------------------------
c     return if matrix is singular (zero pivot) after column interchange
c     -----------------------------------------------------------------
      if ( abs(det) .eq. 0e0 ) go to 600
c     -----------------------------------------------------------------
c     -----------------------------------------------------------------
      pivot(ic) = amax
c     -----------------------------------------------------------------
c     interchange rows to put pivot element on diagonal
c     -----------------------------------------------------------------
      if ( ir .eq. ic ) go to 260
      det = -det
      do 200 k = 1, n
      swap = a(j,k)
      a(j,k) = a(ic,k)
      a(ic,k) = swap
  200 continue
      if ( m .le. 0  ) go to 260
      do 250 k = 1, m
      swap = b(j,k)
      b(j,k) = b(ic,k)
      b(ic,k) = swap
  250 continue
c     -----------------------------------------------------------------
c     divide pivot row by pivot element
c     -----------------------------------------------------------------
  260 do 350 k = 1, n
      if ( k .eq. ic ) a(ic,k) = 1.
      a(ic,k) = a(ic,k)/t
  350 continue
      if ( m .le. 0 ) go to 380
      do 370 k = 1, m
      b(ic,k) = b(ic,k)/t
  370 continue
c     -----------------------------------------------------------------
c     reduce non-pivot rows
c     -----------------------------------------------------------------
  380  do 550 j = 1, n
       if ( j .eq. ic ) go to 550
       t = a(j,ic)
       a(j,ic) = 0e0
       do 450 k = 1, n
       a(j,k) = a(j,k) - a(ic,k)*t
  450 continue
      if ( m .le. 0 ) go to 550
      do 500 k = 1, m
      b(j,k) = b(j,k) - b(ic,k)*t
  500 continue
  550 continue
      i = n
c     -----------------------------------------------------------------
c     interchange colums after all pivot operations have been performed
c     -----------------------------------------------------------------
  600 do 710 il = 1, i
      ii = i + 1 - il
      k = index(ii)/4096
      ic = index(ii) - 4096*k
      if ( k .eq. ic ) go to 710
      do 705 j = 1, n
      swap = a(j,k)
      a(j,k) = a(j,ic)
      a(j,ic) = swap
  705 continue
  710 continue
      return
      end
