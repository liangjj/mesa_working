*deck expad
      function expad ( x, n )
      real *8 x, expad
      IF ( n .eq. 0 ) then
         expad = x
      else
         expad = x*2d0**n
      endif
      return
      end
