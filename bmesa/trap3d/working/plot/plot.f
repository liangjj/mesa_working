      program plot
      implicit integer(a-z)
      common // ia(1)
      real*8 z(1)
      equivalence (ia(1),z(1))
      common /memory/ ioff
c
c
      call drum
      call getscm(need,z,maxcor,'plot',0)
      intoff=wpadti(ioff)
      call mkplot(z(ioff),ia(intoff))
c
c
      stop
      end
