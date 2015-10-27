*deck mkwt.f
c***begin prologue     mkwt
c***date written       930502   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           mkwt, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            calculate atomic weights
c***references         
c
c***routines called
c***end prologue       mkwt
      subroutine mkwt (wtg,wtnc,wtth,wtph,wtleb,wt,scr,nr,nthet,nphi,
     1                 typ,nang,nonsep)
      implicit integer (a-z)
      real*8 wtg, wtnc, wtth, wtph, wt, wtleb, scr
      character*(*) typ
      logical nonsep
      dimension wtg(nr), wtnc(nr,nr-1), wtth(nthet), wtph(nphi), wt(*)
      dimension wtleb(nang), scr(*)
      common /io/ inp, iout
      count=0
      if (.not.nonsep) then
          angcnt=0
          do 10 k=1,nthet
             do 20 l=1,nphi
                angcnt=angcnt+1
                scr(angcnt)=wtth(k)*wtph(l)
   20        continue
   10     continue     
          if (typ.eq.'newton-cotes') then
              do 30 i=1,nr-1
                 do 40 j=1,nr
                    call vscale(wt(count+1),scr,wtnc(j,i),angcnt)
                    count=count+angcnt
   40            continue
   30         continue
          else
              do 50 i=1,nr
                 call vscale(wt(count+1),scr,wtg(i),angcnt)
                 count=count+angcnt
   50         continue
          endif
      else
          if (typ.eq.'newton-cotes') then
              do 60 i=1,nr-1
                 do 70 j=1,nr
                    call vscale(wt(count+1),wtleb,wtnc(j,i),nang)
                    count=count+nang
   70            continue
   60         continue
          else
              do 80 i=1,nr
                 call vscale(wt(count+1),wtleb,wtg(i),nang)
                 count=count+nang
   80         continue
          endif
      endif                                          
      return
      end

