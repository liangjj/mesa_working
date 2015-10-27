*deck scalwt.f
c***begin prologue     scalwt
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           scalwt, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            scale radial weights.
c
c***routines called
c***end prologue       unscal
      subroutine scalwt(r,wtg,wtnc,scr,n,typ)
      implicit integer (a-z)
      real*8 r, wtg, wtnc, scr
      character*(*) typ
      dimension r(n), wtg(n), wtnc(n,n-1), scr(*)
      call vmul(scr,r,r,n)
      if (typ.eq.'newton-cotes') then
          do 10 i=1,n-1
             call vmul(wtnc(1,i),wtnc(1,i),scr,n)
   10     continue
      else
          call vmul(wtg,wtg,scr,n)
      endif                       
      return
      end

