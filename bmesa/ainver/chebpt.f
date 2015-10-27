*deck chebpt
      subroutine chebpt(x,wt,left,right,n,prnt)
c***begin prologue     chebpt
c***date written       940504   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose             
c***                   .
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       chebpt
c
      implicit integer (a-z)
      real*8 x, wt, left, right
      real*8 add, pi, prefac 
      character*80 title
      logical prnt
      dimension x(n), wt(n)
      common /io/ inp, iout
      data pi/3.1415926535897932384d0/
c
c        use trapezoidal points on (0,pi)
c
      add=(right-left)/(n+1)
      prefac = pi/(right-left)
      x(1) = left + add
      wt(1)=add
      do 10 i=2,n
         x(i) = x(i-1) + add
         wt(i) = add
 10   continue
      if(prnt) then
         title='chebshev points'
         call prntfm(title,x,n,1,n,1,iout)
         title='chebshev weights'
         call prntfm(title,wt,n,1,n,1,iout)
      endif
      return
      end















