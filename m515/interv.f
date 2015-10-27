*deck @(#)interv.f	5.1  11/28/95
      subroutine interv(xt,lxt,x,ileft,mflag)
c
c computes largest ileft in (1,lxt) such that xt(ileft) .le. x
      implicit real *8 (a-h,o-z)
      dimension xt(lxt)
      data ilo /1/
c***first executable statement  interv
      ihi = ilo + 1
      if (ihi .lt. lxt)                go to 20
         if (x .ge. xt(lxt))           go to 110
         if (lxt .le. 1)               go to 90
         ilo = lxt - 1
                                       go to 21
   20 if (x .ge. xt(ihi))              go to 40
   21 if (x .ge. xt(ilo))              go to 100
c *** now x .lt. xt(ihi) . find lower bound
   30 istep = 1
   31 ihi = ilo
      ilo = ihi - istep
      if (ilo .le. 1)                  go to 35
      if (x .ge. xt(ilo))              go to 50
      istep = istep*2
                                       go to 31
   35 ilo = 1
      if (x .lt. xt(1))                go to 90
                                       go to 50
c *** now x .ge. xt(ilo) . find upper bound
   40 istep = 1
   41 ilo = ihi
      ihi = ilo + istep
      if (ihi .ge. lxt)                go to 45
      if (x .lt. xt(ihi))              go to 50
      istep = istep*2
                                       go to 41
   45 if (x .ge. xt(lxt))              go to 110
      ihi = lxt
c *** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval
   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)             go to 100
c     note. it is assumed that middle = ilo in case ihi = ilo+1
      if (x .lt. xt(middle))           go to 53
         ilo = middle
                                       go to 50
   53    ihi = middle
                                       go to 50
c *** set output and return
   90 mflag = -1
      ileft = 1
                                       return
  100 mflag = 0
      ileft = ilo
                                       return
  110 mflag = 1
      ileft = lxt
                                       return
      end
