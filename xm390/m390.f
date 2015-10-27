*deck @(#)m390.f	2.1  10/10/91
      program m390
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m390',0)
      intoff=wpadti(ioff)
      call pm390(z(ioff),ia(intoff))
c
c
      stop
      end
