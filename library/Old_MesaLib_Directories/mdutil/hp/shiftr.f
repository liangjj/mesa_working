*deck @(#)shiftr.f	5.1  11/6/94
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
