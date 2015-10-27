*deck cmpre.f
c***begin prologue     cmpre
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
c***end prologue       cmpre
      subroutine cmpre(d,p,q,type,n)
      implicit integer (a-z)
      real*8 d, p, q, fn, fe
      character*(*) type
      dimension d(n), p(n,n), q(n)
      common/io/inp, iout
      if(type.eq.'sine') then
         do 10 i=1,n
            fn = d(i)*p(i,i)
            fe = sin(q(i))
            write(iout,1) q(i), fe, fn
 10      continue   
      elseif(type.eq.'cosine') then
         do 20 i=1,n
            fn = d(i)*p(i,i)
            fe = cos(q(i))
            write(iout,1) q(i), fe, fn
 20      continue   
      elseif(type.eq.'quadratic') then
         do 30 i=1,n
            fn = d(i)*p(i,i)
            fe = q(i)*q(i)
            write(iout,1) q(i), fe, fn
 30      continue   
      elseif(type.eq.'exponential') then
         do 40 i=1,n
            fn = d(i)*p(i,i)
            fe = exp(-q(i))
            write(iout,1) q(i), fe, fn
 40      continue   
      elseif(type.eq.'gaussian') then
         do 50 i=1,n
            fn = d(i)*p(i,i)
            fe = exp(-q(i)*q(i))
            write(iout,1) q(i), fe, fn
 50      continue   
      else
         call lnkerr('error in function type')
      endif
      return
 1    format(/,1x,'x = ',e15.8,1x,'exact = ',e15.8,1x,
     1                            'approximate = ',e15.8)
      end
