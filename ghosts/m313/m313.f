*deck @(#)m313.f	1.1  11/30/90
      program m313
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m313',0)
      intoff=wpadti(ioff)
      call pm313(z(ioff),ia(intoff))
c
c
      stop
      end
