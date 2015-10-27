*deck @(#)recur2.f	5.1  11/6/94
      subroutine recur2(nhi,lalo,lahi,lblo,lbhi,rka,rkb,qq)
      implicit real*8(a-h,o-z)
c
c     controls computation of remaining q(n,la,lb) by recurrence.
c
      common/const/zero,one,two,three,four,five,six,ten
      dimension qq(11,9,9)
      common/io/inp,iout
c
cdir$ ivdep
      fka=one/rka
      fkb=one/rkb
      lahm1=max(lahi-1,1)
      lbhm1=max(lbhi-1,1)
      do 20 lb=lbhi,lblo+2,-1
         lbtru=lb-1
         do 10 n=1,nhi
            qq(n+1,lahi,lb-2)=qq(n+1,lahi,lb)
     $                       +(2*lbtru-1)*fkb*qq(n,lahi,lb-1)
            qq(n+1,lahm1,lb-2)=qq(n+1,lahm1,lb)
     $                        +(2*lbtru-1)*fkb*qq(n,lahm1,lb-1)
   10    continue
   20 continue
c
      do 50 la=lahi,lalo+2,-1
         latru=la-1
         do 40 n=1,nhi
            qq(n+1,la-2,lbhi)=qq(n+1,la,lbhi)
     $                       +(2*latru-1)*fka*qq(n,la-1,lbhi)
            qq(n+1,la-2,lbhm1)=qq(n+1,la,lbhm1)
     $                        +(2*latru-1)*fka*qq(n,la-1,lbhm1)
   40    continue
   50 continue
c
      do 90 lb=lbhi,lblo+2,-1
         lbtru=lb-1
         do 80 la=lahi,lalo+2,-1
            latru=la-1
            do 70 n=1,nhi
               qq(n+1,la-2,lb-2)=qq(n+1,la-2,lb)
     $                          +(2*lbtru-1)*fkb*qq(n,la-2,lb-1)
   70       continue
   80    continue
   90 continue
c
      return
      end
