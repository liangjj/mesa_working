*deck fil3.f
c***begin prologue     fil3
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***references         
c                        
c***routines called    
c***end prologue       fil3
      subroutine fil3(psi0,phix,phiy,phiz,nx,ny,nz,nc)
      implicit integer (a-z)
      real*8 psi0, phix, phiy, phiz
      dimension psi0(nz,ny,nx,nc), phix(nx), phiy(ny), phiz(nz)
      common/io/inp, iout
      do 10 ic=1,nc
         do 20 i=1,nx
            do 30 j=1,ny
               call copy(phiz,psi0(1,j,i,ic),nz)
 30         continue   
 20      continue
         do 40 k=1,nz
            do 50 i=1,nx
               do 60 j=1,ny
                  psi0(k,j,i,ic)= phiy(j)
 60            continue
 50         continue   
 40      continue
         do 70 k=1,nz
            do 80 j=1,ny
               do 90 i=1,nx
                  psi0(k,j,i,ic)=phix(i)
 90            continue   
 80         continue   
 70      continue   
 10   continue   
      return
      end       

