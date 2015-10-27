*deck mkhpsi.f
c***begin prologue     mkhpsi
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            right hand side of inhomogeneous time-dependent
c***                   hamiltonian.
c***                   
c***description        matrix element of hamiltonian times
c***                   initial state.  
c***                   n = nspace = nx one dimension
c***                              = ny*nx two dimensions
c***                              = nz*ny*nx three dimensions
c
c                      note that the vectors are stored exactly as above.
c***references         
c
c***routines called    
c***end prologue       mkhpsi
      subroutine mkhpsi(hx,hy,hz,driver,psi0,vxyz,tmp,n,nx,ny,nz,
     1                  spac,dim,prnt)
      implicit integer (a-z)
      real*8 hx, hy, hz, vxyz, fac
      real*8  driver, psi0, tmp
      logical prnt, spac
      character*80 title
      character*16 fptoc
      dimension hx(nx,nx), hy(ny,ny), hz(nz,nz)
      dimension driver(n), vxyz(n), psi0(n)
      dimension tmp(n,2)
      common/io/inp, iout
c
      call rzero(driver,n)
      if(dim.eq.1) then
         call fildr1(driver,vxyz,psi0,nx)
         call mkh1d(hx,driver,psi0,tmp,nx)
      elseif(dim.eq.2) then
         call fildr2(driver,vxyz,psi0,nx,ny)
         call mkh2d(hx,hy,driver,psi0,tmp,nx,ny)
      elseif(dim.eq.3) then
         call fildr3(driver,vxyz,psi0,nx,ny,nz)        
         call mkh3d(hx,hy,hz,driver,psi0,tmp,nx,ny,nz)
      endif   
      if(prnt) then
         title='driving term'
         call prntrm(title,driver,n,1,n,1,iout)
      endif                     
      return
      end       
