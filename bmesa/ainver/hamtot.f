*deck hamtot.f
c***begin prologue     hamtot
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            time-dependent hamiltonian not including
c***                   non-linear term
c***                   
c***description        total hamiltonian in space and time dimension
c***                   contructed explicitly.  
c***references         
c
c***routines called    
c***end prologue       hamtot
      subroutine hamtot(ham,hx,hy,hz,energy,vxyz,ind,nx,ny,nz,n,dim,prn)
      implicit integer (a-z)
      real*8 hx, hy, hz, vxyz, ham, energy
      logical prn
      character*80 title
      dimension ham(n,n), vxyz(n), ind(*)
      dimension hx(nx,nx), hy(ny,ny), hz(nz,nz)
      common/io/inp, iout      
      call rzero(ham,n*n)
      if(dim.eq.1) then
         call hamfl1(ham,hx,vxyz,n,prn)
      elseif(dim.eq.2) then
         call setd2(ind,ny,nx,n)
         call hamfl2(ham,hx,hy,vxyz,ind,nx,ny,n,prn)
      elseif(dim.eq.3) then
         call setd3(ind,nz,ny,nx,n)
         call hamfl3(ham,hx,hy,hz,vxyz,ind,nx,ny,nz,n,prn)
      else
         call lnkerr('error in dimension')
      endif  
      do 10 i=1,n
         ham(i,i) = ham(i,i) - energy
 10   continue   
      if (prn) then
          title='full space and time-dependent hamiltonian'
          call prntrm(title,ham,n,n,n,n,iout)
      endif     
      return
      end       
