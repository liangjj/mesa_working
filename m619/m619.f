*deck @(#)m619.f	5.1 11/6/94
      program m619
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
      call drum
      call getscm(0,z,maxcor,'m619',0)
      intoff=wpadti(ioff)
      call pm619(z(ioff),ia(intoff))
c
      stop
      end
