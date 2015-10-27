*deck @(#)m105.f	1.1  11/20/92
      program m105
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m105',0)
      intoff=wpadti(ioff)
      call pm105(z(ioff),ia(intoff))
c
c
      stop
      end
