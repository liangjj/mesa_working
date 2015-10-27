*deck @(#)m1990.f	4.1  7/7/93
      program m1990
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
      call drum
      call getscm(need,z,maxcor,'m1990',0)
      intoff=wpadti(ioff)
      call pm1990(z(ioff),ia(intoff),maxcor)
      stop
      end
