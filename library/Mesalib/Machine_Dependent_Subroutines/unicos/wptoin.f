*deck @(#)wptoin.f	5.1  11/6/94
      function wptoin(i)
c
c     number of working precision words in an integer.
c     64 bit integer/64 bit working precision(real*8) version.
c
      implicit integer(a-z)
      integer wptoin
c
c
      wptoin=i
c
c
      return
      end
