*deck %W%  %G%
      program m203
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
      call drum
      call getscm(need,z,maxcor,'m203',0)
      intoff=wpadti(ioff)
      call pm203(z(ioff),ia(intoff))
c
      stop
      end
