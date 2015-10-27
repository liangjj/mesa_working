*deck @(#)recurf.f	5.1  11/6/94
      subroutine recurf(lalo,lahi,lblo,lbhi,rka,rkb)
      implicit real*8(a-h,o-z)
c
      common/ptwtdat/ptpow(50,11),f(50,9,9),pt(50)
      common/ptwt/npts
      dimension fka(50), fkb(50)
c
      parameter (one=1.0d0)
c
c     use the recurrence relations for modified spherical bessel
c     functions to generate the integrand for remaining lamda values.
c
cdir$ ivdep
      lahm1=max(lahi-1,1)
      lbhm1=max(lbhi-1,1)
      do 10 i=1,npts
         fka(i)=one/(rka*pt(i))
         fkb(i)=one/(rkb*pt(i))
   10 continue
c
      do 30 lb=lbhi,lblo+2,-1
         lbtru=lb-1
         do 20 i=1,npts
            f(i,lahi,lb-2)=f(i,lahi,lb)
     $                    +(2*lbtru-1)*fkb(i)*f(i,lahi,lb-1)
            f(i,lahm1,lb-2)=f(i,lahm1,lb)
     $                    +(2*lbtru-1)*fkb(i)*f(i,lahm1,lb-1)
   20    continue
   30 continue
c
      do 60 la=lahi,lalo+2,-1
         latru=la-1
         do 50 i=1,npts
            f(i,la-2,lbhi)=f(i,la,lbhi)
     $                    +(2*latru-1)*fka(i)*f(i,la-1,lbhi)
            f(i,la-2,lbhm1)=f(i,la,lbhm1)
     $                    +(2*latru-1)*fka(i)*f(i,la-1,lbhm1)
   50    continue
   60 continue
c
      do 90 lb=lbhi,lblo+2,-1
         lbtru=lb-1
         do 80 la=lahi,lalo+2,-1
            latru=la-1
            do 70 i=1,npts
               f(i,la-2,lb-2)=f(i,la-2,lb)
     $                       +(2*lbtru-1)*fkb(i)*f(i,la-2,lb-1)
   70       continue
   80    continue
   90 continue
c
c
      return
      end
