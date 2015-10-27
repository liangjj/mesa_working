*deck %W%  %G%
      program m111
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
      call drum
      call getscm(need,z,maxcor,'m111',0)
      intoff=wpadti(ioff)
      call pm111(z(ioff),ia(intoff))
c
      stop
      end
