*deck condens
      subroutine condens ( e, p, sr, si, i, n )
      implicit real *8 (a-h,o-z)
      dimension e(n), p(n), sr(n), si(n)
c
c     eliminates extra energies
c
      pavg = 0.5d0*(p(i) + p(i+1) )
      p(i) = pavg
      sr(i) = 0.5d0*(sr(i) + sr(i+1))
      si(i) = 0.5d0*(si(i) + si(i+1))
      j = i + 2
      do 1 k = j, n
         e(k - 1) = e(k)
         p(k - 1) = p(k)
         sr(k - 1) = sr(k)
         si(k - 1) = si(k)
    1 continue
      n = n - 1
      return
      end
