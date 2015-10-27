*deck nwtprn
      subroutine nwtprn(pt,wt,n)
c***begin prologue     nwtprn
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            print newton-cotes points and weights
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       nwtprt
c
      implicit integer (a-z)
      real*8 pt, wt
      dimension pt(n), wt(n,n-1)
      common /io/ inp, iout
      do 10 i=1,n-1
         write(iout,1) i
         write(iout,2) pt(1), pt(i+1)
         write(iout,3) (wt(j,i),j=1,n)
   10 continue
      return
 1    format(/,1x,'newton-cotes subinterval = ',i2)
 2    format(/,1x,'weights from x(left) = ',e15.8,' to x(right) = ',
     1             e15.8)
 3    format( (/,1x,5e15.8) )
      end















