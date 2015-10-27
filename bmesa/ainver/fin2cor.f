*deck fin2cor.f
c***begin prologue     fin2cor
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
c***end prologue       fin2cor
      subroutine fin2cor(f,c,p,q,type,nf,nc)
      implicit integer (a-z)
      real*8 f, c, p, q, fe
      character*(*) type
      dimension f(nf), c(nf), p(nc,nf), q(nc)
      common/io/inp, iout
      call ebc(c,p,f,nc,nf,1)
      if(type.eq.'sine') then
         do 10 i=1,nc
            fe = sin(q(i))
            write(iout,1) q(i), fe, c(i)
 10      continue   
      elseif(type.eq.'cosine') then
         do 20 i=1,nc
            fe = cos(q(i))
            write(iout,1) q(i), fe, c(i)
 20      continue   
      elseif(type.eq.'quadratic') then
         do 30 i=1,nc
            fe = q(i)*q(i)
            write(iout,1) q(i), fe, c(i)
 30      continue   
      elseif(type.eq.'exponential') then
         do 40 i=1,nc
            fe = exp(-q(i))
            write(iout,1) q(i), fe, c(i)
 40      continue   
      elseif(type.eq.'gaussian') then
         do 50 i=1,nc
            fe = exp(-q(i)*q(i))
            write(iout,1) q(i), fe, c(i)
 50      continue   
      else
         call lnkerr('error in function type')
      endif
      return
 1    format(/,1x,'x = ',e15.8,1x,'exact = ',e15.8,1x,
     1                            'approximate = ',e15.8)
      end
