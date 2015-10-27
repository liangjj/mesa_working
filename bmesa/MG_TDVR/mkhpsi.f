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
c***description        matrix element of negative of ( i d/dt - H ) times
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
      subroutine mkhpsi(hx,hy,hz,driver,psi0,vxyzt,pt,wt,energy,
     1                  n,nx,ny,nz,nt,dim,i0stat,prnt)
      implicit integer (a-z)
      real*8 hx, hy, hz, vxyzt, pt, wt, fac
      real*8  driver, psi0, energy, tmp, sum, diff
      logical prnt, spac
      character*(*) i0stat
      character*80 title
      character*16 fptoc
      dimension hx(nx,nx), hy(ny,ny), hz(nz,nz)
      dimension pt(nt,nt), wt(nt), driver(n,nt,2), vxyzt(n,nt)
      dimension psi0(n,2)
      common/io/inp, iout
      pointer(ptmp,tmp(1))
c
      need=wptoin(2*n)
      call memory(need,ptmp,ngot,'tmp',0)
      call rzero(driver,n*nt*2)
      if(dim.eq.1) then
         call fildr1(driver,vxyzt,psi0,nx,nt)
         call mkh1d(hx,driver,psi0,tmp,nx,nt)
      elseif(dim.eq.2) then
         call fildr2(driver,vxyzt,psi0,nx,ny,nt)
         call mkh2d(hx,hy,driver,psi0,tmp,nx,ny,nt)
      elseif(dim.eq.3) then
         call fildr3(driver,vxyzt,psi0,nx,ny,nz,nt)        
         call mkh3d(hx,hy,hz,driver,psi0,tmp,nx,ny,nz,nt)
      endif   
      do 40 i=1,nt
         fac = wt(i)*pt(i,i)
         call sscal(n,fac,driver(1,i,1),1)
         call sscal(n,fac,driver(1,i,2),1)
 40   continue   
      if(prnt) then
         title='real part of driving term '
         call prntrm(title,driver(1,1,1),n,nt,n,nt,iout)
         title='imaginary part of driving term '
         call prntrm(title,driver(1,1,2),n,nt,n,nt,iout)
      endif               
c      call memory(-ngot,ptmp,idum,'tmp',idum)
      return
 1    format(/,1x,'test sum = ',e15.8,1x,'i = ',i3,1x,'j = ',i3)
      end       
