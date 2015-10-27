*deck fildr2
c***begin prologue     fildr2
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***routines called    
c***end prologue       mkh2d
      subroutine fildr2(driver,vxy,psi0,nx,ny)
      implicit integer (a-z)
      real*8 vxy, driver, psi0
      dimension driver(ny,nx), psi0(ny,nx), vxy(ny,nx)
      common/io/inp, iout
      do 10 i=1,nx
         do 20 j=1,ny
            driver(j,i) = - vxy(j,i)*psi0(j,i)
 20      continue   
 10   continue
      return
      end       
