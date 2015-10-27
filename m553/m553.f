*deck @(#)m553.f	5.1  11/6/94
      program m553
      implicit integer(a-z)
      common // ia(1)
      common /memory/ ioff
      integer wpadti
      real*8 z(1)
      equivalence (ia(1),z(1))
c
c
      call drum
      call getscm(need,z,maxcor,'m553',0)
      intoff=wpadti(ioff)
      call pm553(z(ioff),ia(intoff))
c
c
      stop
      end
