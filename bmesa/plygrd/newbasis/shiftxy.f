*deck shiftxy.f
c***begin prologue     shiftxy
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose             points, weights and polynomials for non (-1,1) shiftxy
c***                   
c***references         
c
c***routines called    
c***end prologue       shiftxy
      subroutine shiftxy(x,xwt,y,ywt,a,b,n)
      implicit integer (a-z)
      real*8 x, xwt, y, ywt, a, b
      dimension x(n), y(n), xwt(n), ywt(n)
      common/io/inp, iout
      call smuls(y,x,a,b,n)
      call smul(ywt,xwt,a,n)
      return
      end       





