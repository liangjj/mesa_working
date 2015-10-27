*deck @(#)shiftr.f	5.1  11/6/94
      function shiftr(i,j)
c
c     shift the word i to the rigth by j bits.
c
      implicit integer(a-z)
      integer shiftr
c
c
      shiftr=rshift(i,j)
c
c
      return
      end
