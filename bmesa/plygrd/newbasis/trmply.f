*deck trmply.f
c***begin prologue     trmply
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
c***end prologue       trmply
      subroutine trmply(pin,pout,ni,nj,mi,mj,which)
      implicit integer (a-z)
      real*8 pin, pout
      character*(*) which
      dimension pin(nj,ni), pout(mj,mi)
      common/io/inp, iout
      if(which.eq.'first') then
         do 10 i=2,ni
            call copy(pin(2,i),pout(1,i-1),mj)
 10      continue   
      elseif(which.eq.'last') then
         do 20 i=1,ni-1
            call copy(pin(1,i),pout(1,i),mj)
 20      continue
      endif   
      return
      end       



