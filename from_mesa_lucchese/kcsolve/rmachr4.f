      subroutine rmachr4(rmach)

      real*4 rmach(5)
c      real*4 two
c      integer nexp, nman

c      nexp = 8
c      nman = 24
c      two = 2

c      rmach(1) = two**(-(2**(nexp-1) -3))
c      rmach(2) = two**(2**(nexp-1) -1)

c      rmach(3) = two**(1-nman)
c      rmach(4) = two**(1-nman)

c      rmach(5) = log10(two)

      rmach(1) = 0.2350989E-37
      rmach(2) = 0.1701412E+39
      rmach(3) = 0.1192093E-06
      rmach(4) = 0.1192093E-06
      rmach(5) = 0.3010300E+00

      return
      end

