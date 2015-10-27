*deck %W%  %G%
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
      subroutine mkwt (wtg,wtnc,wtth,wtph,wtleb,wt,nr,nthet,nphi,typ,
     1                 nang,nonsep)
      implicit integer (a-z)
      real*8 wtg, wtnc, wtth, wtph, wt, wtleb
      character*(*) typ
      logical nonsep
      dimension wtg(nr), wtnc(nr,nr-1), wtth(nthet), wtph(nphi), wt(*)
      dimension wtleb(nang)
      common /io/ inp, iout
      count=0
      if (.not.nonsep) then
          if (typ.eq.'newton-cotes') then
              do 10 i=1,nr-1
                 do 20 j=1,nr
                    do 30 k=1,nthet
                       do 40 l=1,nphi
                          count=count+1
                          wt(count)=wtnc(j,i)*wtth(k)*wtph(l)
   40                  continue
   30               continue
   20            continue
   10         continue
          else
              do 50 i=1,nr
                 do 60 j=1,nthet
                    do 70 l=1,nphi
                       count=count+1
                       wt(count)=wtg(i)*wtth(j)*wtph(l)
   70               continue
   60            continue
   50         continue
          endif
      else
          if (typ.eq.'newton-cotes') then
              do 80 i=1,nr-1
                 do 90 j=1,nr
                    do 100 k=1,nang
                       count=count+1
                       wt(count)=wtnc(j,i)*wtleb(k)
  100               continue
   90            continue
   80         continue
          else
              do 110 i=1,nr
                 do 120 j=1,nang
                    count=count+1
                    wt(count)=wtg(i)*wtleb(j)
  120            continue
  110         continue
          endif
      endif                                          
      return
      end

