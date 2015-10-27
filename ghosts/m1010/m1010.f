*deck @(#)m1010.f	1.1  11/30/90
      program m1010
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
      call drum
      call getscm(need,z,maxcor,'m1010',0)
      intoff=wpadti(ioff)
      call pm1010(z(ioff),ia(intoff))
      stop
      end
