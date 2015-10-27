*deck fildr
c***begin prologue     fildr
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***routines called    
c***end prologue       fildr
      subroutine fildr(driver,v,psi0,n,nt,nc)
      implicit integer (a-z)
      real*8 v, driver, psi0, one
      character*80 title
      dimension driver(n,nt,nc,2), psi0(n,nc,2), v(n,nt,nc,nc)
      common/io/inp, iout
      data one / 1.d0 /  
      do 10 i=1,n
         do 20 j=1,nt
            call sgemv('n',nc,nc,one,v(i,j,1,1),n*nt*nc,
     1                  psi0(i,1,1),n,one,driver(i,j,1,1),n*nt)
            call sgemv('n',nc,nc,one,v(i,j,1,1),n*nt*nc,
     1                  psi0(i,1,2),n,one,driver(i,j,1,2),n*nt)
 20      continue   
 10   continue
c      do 300 i=1,nc
c         title='driver'
c         call prntrm(title,driver(1,1,i,1),n,nt,n,nt,iout)
c	 call prntrm(title,driver(1,1,i,2),n,nt,n,nt,iout)
c 300  continue   
      return
      end       
