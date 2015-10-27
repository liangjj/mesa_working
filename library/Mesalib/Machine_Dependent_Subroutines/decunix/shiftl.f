*deck @(#)shiftl.f	1.1  11/30/90
      function shiftl(i,j)
c
      implicit integer(a-z)
      integer shiftl
c
c
      shiftl=ishft(i,j)
c
      return
      end
