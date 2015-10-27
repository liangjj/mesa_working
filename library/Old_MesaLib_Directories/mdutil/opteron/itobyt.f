*deck @(#)itobyt.f	5.1  11/6/94
      function itobyt(n)
c
c     number of bytes in an integer 
c
      implicit integer(a-z)
      integer itobyt
c
c*************************************
c     64-bit machines (cray)
c     itobyt=8*n
c*************************************
c
c     32-bit integer version
c
      itobyt=4*n
c
c
      return 
      end
