      subroutine unpak ( x, n )
      implicit real*8 ( a - h, o - z )
c      data mask / 077777b /
c
      IF ( x .ne. 0d0 ) THEN
c         n1 = and(mask,shift(x,16))
c         n2 = and(mask,shift(1e0,16))
c         n = n1 - n2 + 1
         n = 1.442695041d0*log(x) + 1.d0
c         unpak = x*2e0**(-n)
      ELSE
         n = 0
c         unpak = 0e0
      END IF
      return
      end
