*deck wnlt2
      logical function wnlt2 (me, mend, ir, factor, tau, scale, wic)
c***begin prologue  wnlt2
c***subsidiary
c***purpose  subsidiary to wnlit
c***library   slatec
c***type      single precision (wnlt2-s, dwnlt2-d)
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
c***see also  wnilt
c***routines called  (none)
c***revision history  (yymmdd)
c   790701  date written
c   890620  code extracted from wnlit and made a subroutine.  (rwc))
c***end prologue  wnlt2
      real             factor, scale(*), tau, wic(*)
      integer ir, me, mend
c
      real             rn, sn, t
      integer j
c
c***first executable statement  wnlt2
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
      wnlt2 = rn .gt. sn*tau**2
      return
      end
