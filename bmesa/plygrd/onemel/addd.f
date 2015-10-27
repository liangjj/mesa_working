*deck @(#)addd.f
      subroutine addd(a,b,c,n)
c***begin prologue     addd
c***date written       
c***revision date      yymmdd  (yymmdd)
c***keywords           diagonal matrix add
c***author             schneider, barry (nsf)
c***source             
c***purpose            vectorized matrix 
c***description
c                      call addd(a,b,c,n)
c                        a       output matrix, (n,n).
c                        b       input matrix, (n,n).
c                        c       input vector, (n).
c
c***references
c***routines called    
c***end prologue       addd
      implicit integer(a-z)
c
      real*8 a(n,n), b(n,n), c(n)
c
      do 10 i=1,n
         a(i,i) = b(i,i) + c(i)      
 10   continue   
c
      return
      end
