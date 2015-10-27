*deck tql1
      subroutine tql1(n,d,e,ierr)
      implicit none
      integer i,j,l,m,n,ii,l1,l2,mml,ierr
      real*8 d(n),e(n)
      real*8 b,c,c2,c3,dl1,el1,f,g,h,p,r,s,s2
      real*8 pythag ,av ,bv
      pythag(av,bv) = sqrt(av**2+bv**2)
      ierr = 0
      if (n .eq. 1) go to 1001
      do 100 i = 2, n
  100 e(i-1) = e(i)
      f = 0.0d0
      b = 0.0d0
      e(n) = 0.0d0
      do 290 l = 1, n
         j = 0
         h = abs(d(l)) + abs(e(l))
         if (b .lt. h) b = h
         do m = l, n
            if (b + abs(e(m)) .eq. b) go to 120
         end do
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + sign(r,p))
         d(l1) = e(l) * (p + sign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
         do 140 i = l2, n
  140    d(i) = d(i) - h
  145    f = f + h
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            if (abs(p) .lt. abs(e(i))) go to 150
            c = e(i) / p
            r = sqrt(c*c+1.0d0)
            e(i+1) = s * p * r
            s = c / r
            c = 1.0d0 / r
            go to 160
  150       c = p / e(i)
            r = sqrt(c*c+1.0d0)
            e(i+1) = s * e(i) * r
            s = 1.0d0 / r
            c = c * s
  160       p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
  200    continue
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         if (b + abs(e(l)) .gt. b) go to 130
  210    p = d(l) + f
         if (l .eq. 1) go to 250
         do ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
         end do
  250    i = 1
  270    d(i) = p
  290 continue
      go to 1001
 1000 ierr = l
 1001 return
      end
