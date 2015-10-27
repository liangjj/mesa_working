*deck @(#)m890.f	1.1  11/30/90
      program m890
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m890',0)
      intoff=wpadti(ioff)
      call pm890(z(ioff),ia(intoff))
c
c
      stop
      end
