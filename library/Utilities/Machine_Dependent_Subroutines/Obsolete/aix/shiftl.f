*deck @(#)shiftl.f	5.1  11/6/94
      function shiftl(i,j)
c
c     shift the word i to the left by j bits.
c
      implicit integer(a-z)
      integer shiftl
c
c
      shiftl=lshift(i,j)
c
      return
      end
