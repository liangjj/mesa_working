*deck trbak1      
      subroutine trbak1(nm,n,a,e,m,z)
      implicit none
      integer i,j,k,l,m,n,nm
      real*8 a(nm,n),e(n),z(nm,m)
      real*8 s
      if (m .eq. 0) go to 200
      if (n .eq. 1) go to 200
      do 140 i = 2, n
         l = i - 1
         if (e(i) .eq. 0.0d0) go to 140
         do 130 j = 1, m
            s = 0.0d0
            do 110 k = 1, l
  110       s = s + a(i,k) * z(k,j)
            s = (s / a(i,l)) / e(i)
            do 120 k = 1, l
  120       z(k,j) = z(k,j) + s * a(i,k)
  130    continue
  140 continue
  200 return
      end
