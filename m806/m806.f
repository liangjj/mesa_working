*deck @(#)pm806.f	5.1  11/6/94
      program m806
c
      implicit integer (a-z)
c
c
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
      call drum
      call getscm(need,z,maxcor,'m806',0)
      intoff=wpadti(ioff)
      call pm806(z(ioff),ia(intoff),maxcor)
c
c
      stop
      end
