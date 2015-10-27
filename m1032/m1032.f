*deck @(#)m1032.f	5.1  11/6/94
      program m1032
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
      call drum
      call getscm(need,z,maxcor,'m1032',0)
      intoff=wpadti(ioff)
      call pm1032(z(ioff),ia(intoff))
      stop
      end
