*deck @(#)m840.f	5.1  11/6/94
      program m840
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m840',0)
      intoff=wpadti(ioff)
      call pm840(z(ioff),ia(intoff))
c
c
      stop
      end