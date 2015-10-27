*deck mkham.f
c***begin prologue     mkham
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***description        
c***                     
c***references         
c
c***routines called    
c***end prologue       mkham
      subroutine mkham(h,h0,v,n)
      implicit integer (a-z)
      real*8 h, h0, v
      dimension h(n,n), h0(n,n), v(n) 
      common/io/inp, iout 
      call copy(h0,h,n*n)
      do 10 i=1,n
         h(i,i) = h0(i,i) + v(i)
 10   continue   
      return
      end       
