*deck fildr1.f
c***begin prologue     fildr1
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***end prologue       fildr1
      subroutine fildr1(driver,vx,psi0,nx)
      implicit integer (a-z)
      real*8 vx, driver, psi0
      character*80 title
      dimension driver(nx), vxt(nx), psi0(nx)
      common/io/inp, iout
      do 10 i=1,nx
         driver(i) = - vxt(i)*psi0(i) 
 10   continue      
      return
      end       
