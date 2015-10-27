*deck @(#)itobyt.f	5.1  11/6/94
      function itobyt(n)
c
c     number of bytes in an integer 
c     64 bit integer version
c
      implicit integer(a-z)
      integer itobyt
c
c
      itobyt=8*n
c
c     ----- 32 bit integer modification -----
c     itobyt=4*n
c     -----
c
c
      return 
      end
