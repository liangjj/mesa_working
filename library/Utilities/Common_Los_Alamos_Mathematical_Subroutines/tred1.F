*deck tred1      
      subroutine tred1(nm,n,a,d,e,e2)
      implicit none
      integer i,j,k,l,n,ii,nm,jp1
      real*8 a(nm,n),d(n),e(n),e2(n)
      real*8 f,g,h,scale
      do 100 i = 1, n
  100 d(i) = a(i,i)
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
         do 120 k = 1, l
  120    scale = scale + abs(a(i,k))
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 290
  140    do k = 1, l
            a(i,k) = a(i,k) / scale
            h = h + a(i,k) * a(i,k)
         end do
         e2(i) = scale * scale * h
         f = a(i,l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         a(i,l) = f - g
         if (l .eq. 1) go to 270
         f = 0.0d0
         do 240 j = 1, l
            g = 0.0d0
            do 180 k = 1, j
  180       g = g + a(j,k) * a(i,k)
            jp1 = j + 1
            if (l .lt. jp1) go to 220
            do 200 k = jp1, l
  200       g = g + a(k,j) * a(i,k)
  220       e(j) = g / h
            f = f + e(j) * a(i,j)
  240    continue
         h = f / (h + h)
         do j = 1, l
            f = a(i,j)
            g = e(j) - h * f
            e(j) = g
            do k = 1, j
               a(j,k) = a(j,k) - f * e(k) - g * a(i,k)
	    end do
         end do
  270    do 280 k = 1, l
  280    a(i,k) = scale * a(i,k)
  290    h = d(i)
         d(i) = a(i,i)
         a(i,i) = h
  300 continue
      return
      end 
