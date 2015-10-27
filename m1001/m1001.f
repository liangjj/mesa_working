*deck @(#)m1001.f	1.1  11/30/90
      program m1001
c
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m1001',0)
      intoff=wpadti(ioff)
      call pm1001(z(ioff),ia(intoff))
c
c
      stop
      end
