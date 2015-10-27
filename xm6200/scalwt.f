*deck %W%  %G%
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
      subroutine scalwt(r,wtg,wtnc,n,typ)
      implicit integer (a-z)
      real*8 r, wtg, wtnc
      character*(*) typ
      dimension r(n), wtg(n), wtnc(n,n-1)
      if (typ.eq.'newton-cotes') then
          do 10 i=1,n-1
             do 20 j=1,n
                wtnc(j,i)=wtnc(j,i)*r(j)*r(j)
   20        continue
   10     continue
      else
          do 30 i=1,n
             wtg(i)=wtg(i)*r(i)*r(i)
   30     continue
      endif                       
      return
      end

