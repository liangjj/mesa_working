      function pnear ( p, px )
      real *8 pnear, p, px, xincr
c
c     p between +1 and -1, add(subtract) multiples of 2
c         until p is closest to px
c
      xincr = sign(2.d0, px-p)
      incr = abs(px-p)*.5d0 + .5d0
      pnear = p + incr*xincr
      return
      end
