*deck @(#)m914.f	2.1  10/10/91
      program m914
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m914',0)
      intoff=wpadti(ioff)
      call izero(ia(intoff),maxcor)
      call pm914(z(ioff),ia(intoff))
c
c
      stop
      end
