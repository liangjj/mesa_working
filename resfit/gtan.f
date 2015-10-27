*deck gtan
      function gtan(c)
      real *8 c, gtan
      gtan=c/(1.d0-abs(c))
      return
      end
