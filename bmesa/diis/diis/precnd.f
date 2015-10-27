*deck precnd.f
c***begin prologue     precnd
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       precnd
      subroutine precnd(a,b,g,n,type)
      implicit integer (a-z)
      real*8 a, b, g, tmp
      character*(*) type
      dimension a(n,n), b(n), g(n)
      common/io/inp, iout
      if(type.eq.'rhs') then
         call copy(b,g,n)
      elseif(type.eq.'inverse-diagonals') then
         do 10 i=1,n
             tmp=1.d0/a(i,i)
             b(i)=b(i)*tmp      
             g(i)=b(i)
             do 20 j=1,n
                a(i,j)=a(i,j)*tmp
 20          continue   
 10      continue   
      endif
      return
      end       
