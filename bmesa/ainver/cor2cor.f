*deck cor2cor.f
c***begin prologue     cor2cor
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polynomials
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***                                                          
c***references         
c
c***routines called    
c***end prologue       cor2cor
      subroutine cor2cor(f,c,p,n)
      implicit integer (a-z)
      real*8 f, c, p
      dimension f(n), c(n), p(n,n)
      common/io/inp, iout
      do 10 i=1,n
         c(i) = f(i)/p(i,i)     
 10   continue   
      return
      end
