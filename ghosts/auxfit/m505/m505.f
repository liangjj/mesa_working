*deck @(#)m505.f	3.1  11/20/92
      program m505
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m505',0)
      intoff=wpadti(ioff)
      call pm505(z(ioff),ia(intoff))
c
c
      stop
      end
