*deck @(#)m604.f	5.1  11/6/94
      program m604
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
      integer wpadti
c
c
      call drum
      call getscm(need,z,maxcor,'m604',0)
      intoff=wpadti(ioff)
      call pm604(z(ioff),ia(intoff),maxcor)
c
c
      stop
      end
