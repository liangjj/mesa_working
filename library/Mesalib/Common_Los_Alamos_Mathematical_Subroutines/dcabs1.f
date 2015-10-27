      function dcabs1(z)
      complex*16 z,zz
      real*8 dcabs1
      real*8 t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = dabs(t(1)) + dabs(t(2))
      return
      end
