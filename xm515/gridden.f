*deck %W% %G%
      subroutine gridden(nbf,nnp,ng,mxgbsiz,d,phi,grad,
     $                   dengrid,dengrad,minesz,scr,phibar,gradbar,
     $                   dmcut,dograd)
c***begin prologue     %M%
c***date written       930521 
c***revision date      %G%
c   february 22, 1995  rlm at lanl
c      revising algorithm to use intermediate sum and grid
c      blocks which allow more accurate estimate of zero contributions.
c   may 13, 1994       rlm at lanl
c      reshaping grad to be (mxgbsiz,nbf,3)
c   february 5, 1994   rlm at lanl
c      extracting a segment of code out of m511/kmatrix
c***keywords           
c***author             russo, tom(lanl) 
c***source             %W% %G%
c***purpose            goes from basis set representation of the density
c                      to grid representation.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       %M%
      implicit none
c     --- input variables -----
      integer nbf,nnp,ng,mxgbsiz,minesz
      logical dograd
      real*8 dmcut
c     --- input arrays (unmodified) ---
      real*8 d(nnp)
      real*8 phi(mxgbsiz,nbf),grad(mxgbsiz,nbf,3)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 dengrid(mxgbsiz)
      real*8 dengrad(mxgbsiz,3)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 scr(minesz,nbf),phibar(nbf),gradbar(nbf)
c     --- local variables ---
      integer kl,mu,nu,coord
      integer inp,iout
      integer ngb,oddblock,ioff,n
      integer gb,i
      real*8 foo,two,zero,sum,test,bfcut
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      logical timeit
c
      parameter (timeit=.false.)
      parameter (zero=0.0d0,two=2.0d+00)
      parameter (bfcut=0.0d0)
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
c     --- divide this atom's grid into smaller blocks.
      ngb=ng/minesz
      oddblock=mod(ng,minesz)
      if(oddblock.ne.0) ngb=ngb+1
      ioff=0
c
c     --- loop over these smaller blocks.
      do 100 gb=1,ngb
         if(gb.ne.ngb.or.oddblock.eq.0) then
            n=minesz
         else
            n=oddblock
         endif
c        --- get mean value of the basis functions on this block.
         do 10 mu=1,nbf
            sum=zero
            do 5 i=1,n
               sum=sum+abs(phi(ioff+i,mu))
    5       continue
            phibar(mu)=sum/float(n)
   10    continue
c        --- loop over basis functions and determine density.
         call rzero(scr,minesz*nbf)
         kl=0
         do 40 mu=1,nbf
c            if(phibar(mu).gt.dmcut) then
            do 30 nu=1,mu
               kl=kl+1
c               kl=(mu*(mu-1)/2)+nu
               foo=d(kl)
               if(mu.ne.nu) foo=foo+foo
c              --- test the contribution of this block.
               test=foo*phibar(mu)*phibar(nu)
c               if (abs(test).gt.bfcut) then
c               if (abs(d(kl)).ge.dmcut) then
               if (abs(test).ge.dmcut) then
                  call saxpy(n,foo,phi(1+ioff,mu),1,scr(1,nu),1)
               endif
 30         continue
c            endif
 40      continue
c
         do 50 nu=1,nbf
c            if(phibar(nu).gt.bfcut) then
            do 45 i=1,n
               dengrid(ioff+i)=dengrid(ioff+i)
     $                        +scr(i,nu)*phi(ioff+i,nu)
   45       continue
c         endif
   50    continue
      
c
c        --- perhaps do the gradient
         if(dograd) then
            do 90 coord=1,3
c              --- get mean value of gradient of basis function
               do 60 mu=1,nbf
                  sum=zero
                  do 55 i=1,n
                     sum=sum+grad(ioff+i,mu,coord)
 55               continue
                  gradbar(mu)=sum/float(n)
 60            continue
               call rzero(scr,minesz*nbf)
               kl=0
               do 80 mu=1,nbf
                  do 70 nu=1,mu
                     kl=kl+1
c                    --- test the contribution of this block.
                     foo=d(kl)+d(kl)
                     test=foo*max(phibar(mu)*gradbar(nu),
     $                    phibar(nu)*gradbar(mu))
c                     test=foo
c                     test=abs(d(kl))
                     if (abs(test).ge.dmcut) then
                        if(mu.eq.nu) then
                           call saxpy(n,foo,grad(1+ioff,mu,coord),1,
     $                                scr(1,nu),1)
                        else
c                          note the fancy footwork. del(phimu) -> scr(,nu)
c                          del(phinu) -> scr(,mu)
                           call saxpy(n,foo,grad(1+ioff,mu,coord),1,
     $                                scr(1,nu),1)
                           call saxpy(n,foo,grad(1+ioff,nu,coord),1,
     $                                scr(1,mu),1)
                        endif
                     endif
   70             continue
   80          continue
c
               do 85 nu=1,nbf
                  do 85 i=1,n
                     dengrad(ioff+i,coord)=dengrad(ioff+i,coord)
     $                                 +scr(i,nu)*phi(ioff+i,nu)
   85          continue
   90       continue
         endif
c        --- increment the offset and do the next block.
         ioff=ioff+n
  100 continue
      if(timeit) then
         call timing(dum4,dum5,dum6)
         write(iout,*) 'time to lay down density and gradients',
     $               dum4-dum1,dum5-dum2,dum6-dum3
      endif
c
c
      return
      end
