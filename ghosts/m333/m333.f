*deck @(#)m333.f	1.1  11/30/90
      program m333
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m333',0)
      intoff=wpadti(ioff)
      call pm333(z(ioff),ia(intoff))
c
c
      stop
      end
