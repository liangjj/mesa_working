      subroutine init(y,bcond,type,neq)
      implicit integer(a-z)
      real*8 y, bcond
      character*80 title
      logical type
      dimension y(neq)
      common/io/inp, iout
c      
      call rzero(y,neq)
      if(type) then
         n=neq/2
         title='initial real and imaginary parts of y'
         if(bcond.ne.0.d0) then
            do 10 i=1,n
               y(i)=1.d0
 10         continue
         endif 
         call prntrm(title,y,n,2,n,2,iout)  
      else
         title='initial y'
         if(bcond.ne.0.d0) then
            do 20 i=1,neq
               y(i)=1.d0
 20         continue
         endif
         call prntrm(title,y,neq,1,neq,1,iout)
      endif         
      return
      end
