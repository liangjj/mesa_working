*deck rabs      
      function rabs(z)
      real*8 rabs, t
      complex*16 z,zz
      dimension t(2)
      equivalence(zz,t)
      zz=z
      rabs=abs(t(1))+abs(t(2))
      return
      end
