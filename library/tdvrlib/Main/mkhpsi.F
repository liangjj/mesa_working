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
      subroutine mkhpsi(h1,h2,h3,driver,psi0,v,pt,wt,energy,
     1                  n,nc,nd,nt,dim,prn)
      implicit integer (a-z)
      real*8 h1, h2, h3, v, pt, wt, fac
      real*8  driver, psi0, energy, tmp, sum, diff
      logical prn, spac
      character*80 title
      character*16 fptoc
      character*2 itoc
      dimension nd(dim)
      dimension h1(nd(1),nd(1)), h2(nd(2),nd(2)), h3(nd(3),nd(3))
      dimension pt(nt,nt), wt(nt), driver(n,nt,nc,2)
      dimension v(n,nt,nc,nc), psi0(n,nc,2)
      common/io/inp, iout
      pointer(ptmp,tmp(1))
c
      need=wptoin(2*n*nc)
      call getmem(need,ptmp,ngot,'tmp',0)
      call rzero(driver,n*nt*2*nc)
      call fildr(driver,v,psi0,n,nt,nc)
c
      if(dim.eq.1) then
         call mkh1d(h1,driver,psi0,tmp,nd(1),nt,2*nc)
      elseif(dim.eq.2) then
         call mkh2d(h1,h2,driver,psi0,tmp,nd,nt,2*nc)
      elseif(dim.eq.3) then
         call mkh3d(h1,h2,h3,driver,psi0,tmp,nd,nt,2*nc)
      endif   
c      do 60 i=1,nc
c         title='after mkh1d'
c         call prntrm(title,driver(1,1,i,1),n,nt,n,nt,iout)
c         call prntrm(title,driver(1,1,i,2),n,nt,n,nt,iout)
c 60   continue   
      do 40 i=1,nt
         fac = wt(i)*pt(i,i)
         do 50 j=1,nc
            call sscal(n,fac,driver(1,i,j,1),1)
            call sscal(n,fac,driver(1,i,j,2),1)
 50      continue   
 40   continue   
      if(prn) then
         do 100 ic=1,nc
            title='real part of driving term channel = '//itoc(ic)
            call prntrm(title,driver(1,1,ic,1),n,nt,n,nt,iout)
            title='imaginary part of driving term channel = '//itoc(ic)
            call prntrm(title,driver(1,1,ic,2),n,nt,n,nt,iout)
 100     continue   
      endif               
c      call getmem(-ngot,ptmp,idum,'tmp',idum)
      return
 1    format(/,1x,'test sum = ',e15.8,1x,'i = ',i3,1x,'j = ',i3)
      end       
