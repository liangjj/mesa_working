*deck @(#)moden.f	5.1  11/28/95
      subroutine moden(nbf,nmo,ng,mxgbsiz,c,phi,grad,
     $                 dengrid,dengrad,minesz,scr,scr2,
     $                 phibar,gradbar,bfcut,dograd)
c***begin prologue     moden.f
c***date written       950514 
c***revision date      11/28/95
c   may 14, 1995       rlm at lanl
c      rewriting gridden to form molecular orbitals on the grid,
c      then square them to get the density. 
c   february 22, 1995  rlm at lanl
c      revising algorithm to use intermediate sum and grid
c      blocks which allow more accurate estimate of zero contributions.
c   may 13, 1994       rlm at lanl
c      reshaping grad to be (mxgbsiz,nbf,3)
c   february 5, 1994   rlm at lanl
c      extracting a segment of code out of m511/kmatrix
c***keywords           
c***author             russo, tom(lanl) 
c***source             @(#)moden.f	5.1 11/28/95
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
c***end prologue       moden.f
      implicit none
c     --- input variables -----
      integer nbf,nmo,ng,mxgbsiz,minesz
      logical dograd
      real*8 bfcut
c     --- input arrays (unmodified) ---
      real*8 c(nbf,nmo)
      real*8 phi(mxgbsiz,nbf),grad(mxgbsiz,nbf,3)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 dengrid(mxgbsiz)
      real*8 dengrad(mxgbsiz,3)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 scr(minesz,nmo),scr2(minesz,nmo),phibar(nbf),gradbar(nbf)
c     --- local variables ---
      integer mu,coord,mo
      integer inp,iout
      integer ngb,oddblock,ioff,n
      integer gb,i
      real*8 two,zero,sum,test
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
      logical timeit
c
      parameter (timeit=.false.)
      parameter (zero=0.0d0,two=2.0d+00)
c
      common/io/inp,iout
c
c     --- form density and gradients on the grid for this atom,
c         also the gradient invariants if requested.
c
c         dengrid(i)=sum over occupied mo: S(mu)[c(mu,mo)*phi(i,mu)]**2
c         dengrad(i)=sum over occupied mo: 2*psi(mo)*gradpsi(mo)
c                                   = 2*psi(mo)*S(mu)[c(mu,mo)*grad(phi(i,mu)]
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
c        --- get mean absolute value of the basis functions on this block.
         do 10 mu=1,nbf
            sum=zero
            do 5 i=1,n
               sum=sum+abs(phi(ioff+i,mu))
    5       continue
            phibar(mu)=sum/float(n)
   10    continue
c        --- prepare the mo on the grid.
         call rzero(scr,minesz*nmo)
         do 40 mu=1,nbf
             if(phibar(mu).ge.bfcut) then
               do 30 mo=1,nmo
c                 --- test the contribution of this block.
                  test=c(mu,mo)*phibar(mu)
                  if (abs(test).ge.bfcut) then
                     call saxpy(n,c(mu,mo),phi(1+ioff,mu),1,scr(1,mo),1)
                  endif
   30          continue
             endif
   40    continue
c
c        --- square and accumulate to get density.
         do 50 mo=1,nmo
            do 45 i=1,n
               dengrid(ioff+i)=dengrid(ioff+i)+scr(i,mo)*scr(i,mo)
   45       continue
   50    continue
      
c
c        --- perhaps do the gradient
         if(dograd) then
            do 90 coord=1,3
c              --- get mean absolute value of gradient of basis function
               do 60 mu=1,nbf
                  sum=zero
                  do 55 i=1,n
                     sum=sum+abs(grad(ioff+i,mu,coord))
   55             continue
                  gradbar(mu)=sum/float(n)
   60          continue
               call rzero(scr2,minesz*nmo)
               do 80 mu=1,nbf
                  if(gradbar(mu).ge.bfcut) then
                     do 70 mo=1,nmo
c                       --- test the contribution of this block.
                        test=c(mu,mo)*gradbar(mu)
                        if (abs(test).ge.bfcut) then
                           call saxpy(n,c(mu,mo),grad(1+ioff,mu,coord),
     $                                1,scr2(1,mo),1)
                        endif
   70                continue
                  endif
   80          continue
c
               do 85 mo=1,nmo
                  do 85 i=1,n
                     dengrad(ioff+i,coord)=dengrad(ioff+i,coord)
     $                                 +two*scr(i,mo)*scr2(i,mo)
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
