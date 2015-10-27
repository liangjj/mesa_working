*deck @(#)m403.f	5.1 11/6/94
      program m403
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcorr,'m403',0)
      intoff=wpadti(ioff)
      call pm403(z(ioff),ia(intoff))
      call chainx(0)
c
c
      stop
      end
