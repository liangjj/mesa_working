*deck @(#)wptbyt.f	1.1  9/6/91
      function wptbyt(n)
c
c number of bytes in a working-precision number (real*8)
c
      implicit integer(a-z)
      integer wptbyt
c
c
      wptbyt=8*n
c
c
      return
      end
