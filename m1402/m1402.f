*deck @(#)m1402.f	5.1  11/6/94
      program m1402
      common // ia(1)
      dimension z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
      integer wpadti
c
c
      call drum
      call getscm(need,z,maxcor,'m1402',0)
      intoff=wpadti(ioff)
      call pm1402(z(ioff),ia(intoff),maxcor)
      call chainx(0)
c
c
      stop
      end
