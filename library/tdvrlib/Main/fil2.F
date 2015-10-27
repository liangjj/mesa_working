*deck fil2.f
c***begin prologue     fil2
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***references         
c                        
c***routines called    
c***end prologue       fil2
      subroutine fil2(psi0,phix,phiy,nx,ny,nc)
      implicit integer (a-z)
      real*8 psi0, phix, phiy
      dimension psi0(ny,nx,nc), phix(nx), phiy(ny)
      common/io/inp, iout
      do 10 ic=1,nc
         do 20 i=1,nx
            call copy(phiy,psi0(1,i,ic),ny)
 20      continue
         do 30 i=1,ny
            do 40 j=1,nx
               psi0(i,j,ic)= phix(j)
 40         continue
 30      continue   
 10   continue   
      return
      end       

