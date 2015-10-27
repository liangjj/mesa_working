*deck @(#)m812.f	5.1  11/6/94
      program m812
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m812',0)
      intoff=wpadti(ioff)
      call pm812(z(ioff),ia(intoff))
c
c
      stop
      end
