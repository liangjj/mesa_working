*deck @(#)gridrho.f	5.1 11/6/94
      subroutine gridrho(nbf,nnp,ng,mxgbsiz,d,phi,grad,
     $                  dengrid,dengrad,scr,dmcut,dograd)
c***begin prologue     gridrho.f
c***date written       940429 
c***revision date      11/6/94      
c***keywords           
c***author             martin, richard(lanl) 
c***source             @(#)gridrho.f	5.1   11/6/94
c***purpose            goes from basis set representation of the density
c                      to grid representation.
c***description
c    
c
c***references
c
c***routines called
c
c***end prologue       gridrho.f
      implicit none
c     --- input variables -----
      integer nbf,ng,mxgbsiz
      logical dograd
c     --- input arrays (unmodified) ---
      real*8 c(nbf,nmo)
      real*8 phi(mxgbsiz,nbf),grad(mxgbsiz,nbf,3)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 dengrid(mxgbsiz)
      real*8 dengrad(mxgbsiz,3)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 scr(mxgbsiz)
c     --- local variables ---
      integer inp,iout
      real*8 zero,one,two
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      logical timeit
c
      parameter (timeit=.false.)
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00)
c
      common/io/inp,iout
c
c     --- form density and gradients on the grid for this atom,
c         also the gradient invariants if requested.
c
c         dengrid(i)=sum over mu,nu P(mu,nu)*phi(i,mu)*phi(i,nu)
c         dengrad(i)=sum over mu,nu P(mu,nu)
c                                  *[ grad(phi(i,mu))*phi(i,nu)
c                                    +phi(i,mu)*grad(phi(i,nu))]
      if(timeit) then
         call timing(dum1,dum2,dum3)
      endif
c
c     --- form the molecular orbitals on the grid.
c         psi(grid,i)=sum(mu) phi(grid,mu)*c(mu,i)
      call sgemm('n','n',ng,nocc,nbf,one,phi,mxgbsiz,
     $            c,nbf,zero,psi,ng)
c     --- square to get the density.
      call vmul(scr,psi,psi,ng*nocc)
c     --- include the occupation number, and accumulate.
c         first the doubly occupied orbitals.
c     are these really included here?
      do 10 i=1,nbe
         call vwxs(dengrid,dengrid,scr(1,i),two,+1,ng)      
   10 continue
c     and the open shells.
      do 20 i=nbe+1,nae
         call vadd(dengrid,dengrid,scr(1,i),ng)
   20 continue
c
c     --- now the gradient of the density
      if (dograd) then
         do 100 coord=1,3
            call sgemm('n','n',ng,nocc,nbf,one,phi,mxgbsiz,
     $            c,nbf,zero,psi,ng)
c     --- square to get the density.
      call vmul(scr,psi,psi,ng*nocc)
c     --- include the occupation number, and accumulate.
c         first the doubly occupied orbitals.
c     are these really included here?
      do 10 i=1,nbe
         call vwxs(dengrid,dengrid,scr(1,i),two,+1,ng)      
   10 continue
c     and the open shells.
      do 20 i=nbe+1,nae
         call vadd(dengrid,dengrid,scr(1,i),ng)
   20 continue





         do 10 coord=1,3
c    would like to have grdient in different form.
            call ebc(scr,grad(1,coord,nu),ng)
            call vwxs(dengrad(1,coord),dengrad(1,coord),
     $                scr,foo,+1,ng)
   10    continue
      endif
c
c
      kl=0
      do 40 mu=1,nbf
         do 30 nu=1,mu
            kl=kl+1
            if (abs(d(kl)).lt.dmcut) goto 30
               if (mu.eq.nu) then
                  call vmul(scr,phi(1,mu),phi(1,nu),ng)
                  call vwxs(dengrid,dengrid,scr,d(kl),+1,
     $                 ng)
                  foo=two*d(kl)
                  if (dograd) then
                     do 10 coord=1,3
                        call vmul(scr,phi(1,mu),grad(1,coord,nu),ng)
                        call vwxs(dengrad(1,coord),dengrad(1,coord),
     $                            scr,foo,+1,ng)
   10                continue
                  endif
               else
                  foo=two*d(kl)
                  call vmul(scr,phi(1,mu),phi(1,nu),ng)
                  call vwxs(dengrid,dengrid,scr,foo,+1,ng)
                  if (dograd) then
                     do 20 coord=1,3
                        call vmul(scr,phi(1,mu),grad(1,coord,nu),ng)
                        call vwxy(scr,scr,phi(1,nu),grad(1,coord,mu),
     $                            +1,ng)
                        call vwxs(dengrad(1,coord),dengrad(1,coord),
     $                            scr,foo,+1,ng)
   20                continue
                  endif
               endif
 30         continue 
 40      continue 
      if(timeit) then
         call timing(dum4,dum5,dum6)
         write(iout,*) 'time to lay down density and gradients',
     $               dum4-dum1,dum5-dum2,dum6-dum3
      endif
c
c
      return
      end
