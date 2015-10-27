*deck vec2di.f
      subroutine vec2di(a,b,c,oper,n)
c***begin prologue     vec2di
c***date written       
c***revision date      yymmdd  (yymmdd)
c***keywords           diagonal matrix add
c***author             schneider, barry (nsf)
c***source             
c***purpose            vectorized matrix 
c***description
c                      call vec2di(a,b,c,oper,n)
c                        a       output matrix, (n,n).
c                        b       input matrix, (n,n).
c                        c       input vector, (n).
c
c***references
c***routines called    
c***end prologue       vec2di
      implicit integer(a-z)
c
      real*8 a(n,n), b(n,n), c(n)
      character*(*) oper
      common/io/inp, iout
c
      if(oper.eq.'add') then
         do 10 i=1,n
            a(i,i) = b(i,i) + c(i)      
 10      continue
      elseif(oper.eq.'subtract') then
         do 20 i=1,n
            a(i,i) = b(i,i) - c(i)      
 20      continue
      elseif(oper.eq.'multiply') then
         do 30 i=1,n
            a(i,i) = b(i,i)*c(i)      
 30      continue
      elseif(oper.eq.'divide') then
         do 40 i=1,n
            a(i,i) = b(i,i)/c(i)      
 40      continue
      else
         call lnkerr('bad operation in vec2di')
      endif	                                             
c
      return
      end
