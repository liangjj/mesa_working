*deck @(#)m551.f	1.2  5/30/91
      program m551
      implicit integer (a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'m551',0)
      intoff=wpadti(ioff)
      call pm551(z(ioff),ia(intoff))
c
c
      stop
      end
