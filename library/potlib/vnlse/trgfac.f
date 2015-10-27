*deck trgfac.f
c***begin prologue     trgfac
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           copy
c***author             schneider, barry (nsf)
c***source             
c***purpose            trig factors
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       trgfac
      subroutine trgfac(vout,scale,omega,q,n,type)
      implicit integer (a-z)
      real*8 vout, scale, omega, q
      character*(*) type
      dimension vout(n), q(n)
      common/io/inp, iout
      if(type.eq.'cosine') then
         do 10 i=1,n
            vout(i) = scale*cos(omega*q(i))
 10      continue
      elseif(type.eq.'sine') then
         do 20 i=1,n
            vout(i) = scale*sin(omega*q(i))
 20      continue
      else
         call lnkerr('error in trig factor')
      endif
      return   
      end       


