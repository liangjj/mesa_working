      function expad ( x, n )
      implicit real*8 ( a - h, o - z )
      if ( n .eq. 0 ) then
         expad = x
      else
         expad = x*2d0**n
      endif
      return
      end
