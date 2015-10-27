*deck fil1.f
c***begin prologue     fil1
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***references         
c                        
c***routines called    
c***end prologue       fil1
      subroutine fil1(psi0,phi,n,nc)
      implicit integer (a-z)
      real*8 psi0, phi
      dimension psi0(n,nc), phi(n)
      common/io/inp, iout
      do 10 ic=1,nc
         call copy(phi,psi0(1,ic),n)
 10   continue   
      return
      end       

