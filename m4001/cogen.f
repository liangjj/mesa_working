*deck @(#)cogen.f	5.1 11/6/94
      subroutine cogen
      implicit double precision(a-h,o-z)
      common zah(31),den(31)
      dimension cfr(30),zp(31),dp(31)
      data zp/31*0.0d0/,dp/31*0.0d0/
      cfr(1)=1.0d0
      do 10 i=2,30
      im1=i-1
   10 cfr(i)=2.0d0/(cfr(im1)*dfloat(im1))
      do 30 i=1,31
      zah(i)=0.0d0
   30 den(i)=0.0d0
      zp(1)=1.0d0
      den(1)=1.0d0
      do 100 n=1,30
      np1=n+1
      do 90 is=1,n
      i=n-is+2
      im1=i-1
      a=zah(i)
      zah(i)=cfr(n)*zah(im1)+zp(i)
      zp(i)=a
      a=den(i)
      den(i)=cfr(n)*den(im1)+dp(i)
   90 dp(i)=a
      a=zah(1)
      zah(1)=zp(1)
      zp(1)=a
      a=den(1)
      den(1)=dp(1)
  100 dp(1)=a
      do 150 i=1,31
      zp(i)=zah(i)
  150 dp(i)=den(i)
      do 200 i=1,31
      j=32-i
      zah(j)=zp(i)
  200 den(j)=dp(i)
      return
      end
