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
      subroutine convt(x,wts,lower,upper,norm0,n,type,prnwpt)
c
      implicit integer (a-z)
      real*8 x, wts, lower, upper, norm0
      character*(*) type
      character*80 title
      logical prnwpt
      common/io/inp, iout 
      dimension x(n), wts(n)
      if (type.ne.'hermite'.or.type.ne.'laguerre') then
          norm0=upper-lower
          do 10 i=1,n
             x(i) = .5d0*( (upper-lower)*x(i) + (upper+lower) )
             wts(i) =.5d0*wts(i)*( upper- lower )
 10       continue
      endif 
      if(prnwpt) then
         title='actual points'
         call prntrm(title,x,n,1,n,1,iout)
         title='actual weights'
         call prntrm(title,wts,n,1,n,1,iout)
      endif
      return
      end       
