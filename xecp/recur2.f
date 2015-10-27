*deck %W%  %G%
      subroutine recur2(n,lalo,lahi,lblo,lbhi,rka,rkb,qq)
      implicit integer(a-z)
c
c     controls computation of remaining q(n,la,lb) by recurrence.
c     this recursion for the type 2 integral is essentially the
c     downward recursion for the modified spherical bessel function.
c        m(l-2,x)=m(l,x)+(2l-1)*m(l-1,x)/x
c
c     ----- arguments unchanged -----
      integer n,lalo,lahi,lblo,lbhi
      real*8 rka,rkb
c     ----- arguments returned -----
      real*8 qq(0:6,0:6,0:6)
c     ----- local variables -----
      real*8 fka,fkb
      real*8 one
      parameter (one=1.0d+00)
c
cdir$ ivdep
      fka=one/rka
      fkb=one/rkb
      lahm1=max(lahi-1,0)
      lbhm1=max(lbhi-1,0)
c
c     ----- recur down on lambda bar for starting values of lambda and all n
      do 20 lb=lbhi,lblo+2,-1
         do 10 nabc=1,n
            qq(nabc,lahi,lb-2)=qq(nabc,lahi,lb)
     $                       +(2*lb-1)*fkb*qq(nabc-1,lahi,lb-1)
            qq(nabc,lahm1,lb-2)=qq(nabc+1,lahm1,lb)
     $                        +(2*lb-1)*fkb*qq(nabc-1,lahm1,lb-1)
   10    continue
   20 continue
c
c     ----- recur down on lambda for starting values of lambda bar and all n 
      do 50 la=lahi,lalo+2,-1
         do 40 nabc=1,n
            qq(nabc,la-2,lbhi)=qq(nabc,la,lbhi)
     $                       +(2*la-1)*fka*qq(nabc-1,la-1,lbhi)
            qq(nabc,la-2,lbhm1)=qq(nabc,la,lbhm1)
     $                        +(2*la-1)*fka*qq(nabc-1,la-1,lbhm1)
   40    continue
   50 continue
c
c     ----- recur down on lambda and lambda bar.
      do 90 lb=lbhi,lblo+2,-1
         do 80 la=lahi,lalo+2,-1
            do 70 nabc=1,n
               qq(nabc,la-2,lb-2)=qq(nabc,la-2,lb)
     $                          +(2*lb-1)*fkb*qq(nabc-1,la-2,lb-1)
   70       continue
   80    continue
   90 continue
c
      return
      end
