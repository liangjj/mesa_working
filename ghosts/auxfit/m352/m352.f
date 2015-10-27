*deck @(#)m352.f	1.2  4/28/92
      program m352
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m352',0)
      intoff=wpadti(ioff)
      call pm352(z(ioff),ia(intoff))
      call chainx(0)
c
c
      stop
      end
