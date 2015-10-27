*deck @(#)m916.f	2.1  10/10/91
      program m916
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m916',0)
      intoff=wpadti(ioff)
      call izero(ia(intoff),maxcor)
      call pm916(z(ioff),ia(intoff))
c
c
      stop
      end
