*deck @(#)m824.f	5.1  11/6/94
      program m824
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m824',0)
      intoff=wpadti(ioff)
      call pm824(z(ioff),ia(intoff))
c
c
      stop
      end
