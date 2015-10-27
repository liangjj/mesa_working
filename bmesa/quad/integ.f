*deck integ
      subroutine integ(f,wt,intgl,n)
c***begin prologue     integ
c***date written       930608   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            integrate a function
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       integ
c
      implicit integer (a-z)
      real*8 f, wt, intgl
      dimension f(n), wt(n,n-1), intgl(0:n-1)
      common /io/ inp, iout
      call rzero(intgl(0),n)
      do 10 i=1,n-1
         intgl(i)=intgl(i-1)
         do 20 j=1,n
            intgl(i)=intgl(i)+wt(j,i)*f(j)
   20    continue
   10 continue
      return
      end















