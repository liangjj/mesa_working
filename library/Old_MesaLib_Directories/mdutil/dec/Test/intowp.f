*deck @(#)intowp.f	5.1  11/6/94
      function intowp(i)
      implicit integer(a-z)
      integer intowp
c
c***************************************
c     for 64 bit machines (cray)
c     intowp=i
c***************************************
c     version for 32 bit machines.
c
      intowp=2*i
c
c
      return
      end
