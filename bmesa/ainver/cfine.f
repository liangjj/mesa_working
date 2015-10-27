*deck cfine.f.f
c***begin prologue     cfine.f
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
c***end prologue       cfine.f
      subroutine cfine(f,p,q,wt,type,n)
      implicit integer (a-z)
      real*8 f, p, q, wt
      character*(*) type
      dimension f(n), p(n,n), q(n), wt(n)
      common/io/inp, iout
      if(type.eq.'sine') then
         do 10 i=1,n
            f(i) = sin(q(i))*wt(i)*p(i,i)
 10      continue   
      elseif(type.eq.'cosine') then
         do 20 i=1,n
            f(i) = cos(q(i))*wt(i)*p(i,i)
 20      continue   
      elseif(type.eq.'quadratic') then
         do 30 i=1,n
            f(i) = q(i)*q(i)*wt(i)*p(i,i)
 30      continue   
      elseif(type.eq.'exponential') then
         do 40 i=1,n
            f(i) = exp(-q(i))*wt(i)*p(i,i)
 40      continue   
      elseif(type.eq.'gaussian') then
         do 50 i=1,n
            f(i) = exp(-q(i)*q(i))*wt(i)*p(i,i)
 50      continue   
      else
         call lnkerr('error in function type')
      endif
      return
      end
