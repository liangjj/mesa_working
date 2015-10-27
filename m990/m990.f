*deck @(#)m990.f	1.1  11/30/90
      program m990
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
      integer wpadti
c
c
      call drum
      call getscm(need,z,maxcor,'m990',0)
      intoff=wpadti(ioff)
      call pm990(z(ioff),ia(intoff))
c
c
      stop
      end
