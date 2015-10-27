*deck convt.f
c***begin prologue     convt
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                   
c
c***references         
c
c***routines called    
c***end prologue       convt
      subroutine convt(x,wt,lower,upper,norm0,n,prnwpt)
c
      implicit integer (a-z)
      real*8 x, wt, lower, upper, norm0
      character*80 title
      logical prnwpt
      common/io/inp, iout 
      dimension x(n), wt(n)
      norm0=upper-lower
      do 10 i=1,n
         x(i) = .5d0*( (upper-lower)*x(i) + (upper+lower) )
         wt(i) =.5d0*wt(i)*( upper- lower )
 10   continue
      if(prnwpt) then
         title='actual points'
         call prntfm(title,x,n,1,n,1,iout)
         title='actual weights'
         call prntfm(title,wt,n,1,n,1,iout)
      endif
      return
      end       
