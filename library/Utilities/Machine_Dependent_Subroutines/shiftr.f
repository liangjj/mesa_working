*deck @(#)shiftr.f	1.1  11/30/90
      function shiftr(i,j)
c
      implicit integer(a-z)
      integer shiftr
c
c
      shiftr=ishft(i,-j)
c
c
      return
      end
