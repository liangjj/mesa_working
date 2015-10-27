*deck dwnlt2
      logical function dwnlt2 (me, mend, ir, factor, tau, scale, wic)
c***begin prologue  dwnlt2
c***subsidiary
c***purpose  subsidiary to wnlit
c***library   slatec
c***type      double precision (wnlt2-s, dwnlt2-d)
c***author  hanson, r. j., (snla)
c           haskell, k. h., (snla)
c***description
c
c     to test independence of incoming column.
c
c     test the column ic to determine if it is linearly independent
c     of the columns already in the basis.  in the initial tri. step,
c     we usually want the heavy weight alamda to be included in the
c     test for independence.  in this case, the value of factor will
c     have been set to 1.e0 before this procedure is invoked.
c     in the potentially rank deficient problem, the value of factor
c     will have been set to alsq=alamda**2 to remove the effect of the
c     heavy weight from the test for independence.
c
c     write new column as partitioned vector
c           (a1)  number of components in solution so far = niv
c           (a2)  m-niv components
c     and compute  sn = inverse weighted length of a1
c                  rn = inverse weighted length of a2
c     call the column independent when rn .gt. tau*sn
c
c***see also  dwnlit
c***routines called  (none)
c***revision history  (yymmdd)
c   790701  date written
c   890620  code extracted from wnlit and made a subroutine.  (rwc))
c   900604  dp version created from sp version.  (rwc)
c***end prologue  dwnlt2
      double precision factor, scale(*), tau, wic(*)
      integer ir, me, mend
c
      double precision rn, sn, t
      integer j
c
c***first executable statement  dwnlt2
      sn = 0.e0
      rn = 0.e0
      do 10 j=1,mend
         t = scale(j)
         if (j.le.me) t = t/factor
         t = t*wic(j)**2
c
         if (j.lt.ir) then
            sn = sn + t
         else
            rn = rn + t
         endif
   10 continue
      dwnlt2 = rn .gt. sn*tau**2
      return
      end
