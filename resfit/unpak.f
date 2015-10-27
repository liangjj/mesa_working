      function unpak ( x, n )
      data mask / 077777b /
c
      IF ( x .ne. 0e0 ) THEN
         n1 = and(mask,shift(x,16))
         n2 = and(mask,shift(1e0,16))
         n = n1 - n2 + 1
         unpak = x*2e0**(-n)
      ELSE
         n = 0
         unpak = 0e0
      END IF
      return
      end
