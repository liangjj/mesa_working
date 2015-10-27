*deck @(#)fmkopen.f	5.1  11/28/95
      subroutine fmkopen(ngtmp,mxgbsiz,minesz,nbf,nnp,becke,lyp,
     $                  tmpgwt,scr,scr2,
     $                  queue,tea,fout,phi,grad,phibar,gradbar,
     $                  dengrada,dengradb,
     $                  kmat,kmcut)
c***begin prologue     fmkopen.f
c***date written       950224  
c***revision date      4/18/95      
c
c***keywords           
c***author             russo, t.v. (lanl)
c***source             @(#)fmkopen.f	5.1   11/28/95
c***purpose            forms the closed shell exchange matrix.   
c***description     
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       fmkopen.f
      implicit none
c     --- input variables -----
      integer ngtmp,nbf,nnp,mxgbsiz,minesz
      logical becke,lyp
      real*8 kmcut
c     --- input arrays (unmodified) ---
      real*8 fout(mxgbsiz,6),phi(mxgbsiz,nbf)
      real*8 grad(mxgbsiz,nbf,3),dengrada(mxgbsiz,3),dengradb(mxgbsiz,3)
      real*8 tmpgwt(ngtmp)
c     --- input arrays (scratch) ---
      real*8 scr(minesz),scr2(minesz)
      real*8 queue(minesz,3),tea(minesz,nbf)
      real*8 phibar(nbf),gradbar(nbf,3)
c     --- output arrays ---
      real*8 kmat(nnp,2)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer mu,nu,i
      integer ngb,oddblock,ioff,gb,n,kl,coord,exch
      integer inp,iout
      real*8 sumo,opbar,sum,test,testnu,test2,test3,test4,test5
      real*8 sdot,pmax,gmax(3)
      real*8 zero
c
      parameter (zero=0.0d0)
      common/io/inp,iout
c
c     --- using the pieces of the functional in fout, form the exchange matrix.
c        fout(*,2) = df/d(rhoa)
c        fout(*,3) = df/d(gaa)
c        fout(*,4) = df/d(rhob)
c        fout(*,5) = df/d(gbb)
c        fout(*,6) = df/d(gab)
c
c     --- divide the current grid block into even smaller blocks.
      ngb=ngtmp/minesz
      oddblock=mod(ngtmp,minesz)
      if(oddblock.ne.0) ngb=ngb+1
      ioff=0
      kmcut=zero
c
c     --- loop over these smaller blocks.
      do 100 gb=1,ngb
         if(gb.ne.ngb.or.oddblock.eq.0) then
            n=minesz
         else
            n=oddblock
         endif
c
c        --- get mean value of the basis functions on this block.
         pmax=zero
         do 20 mu=1,nbf
            sum=zero
            do 15 i=1,n
               sum=sum+abs(phi(ioff+i,mu))
   15       continue
            phibar(mu)=sum/float(n)
            pmax=max(pmax,phibar(mu))
   20    continue
c
c        --- mean value of grad(bf) on this block
         if(becke.or.lyp) then
            do 25 coord=1,3
               gmax(coord)=zero
               do 24 mu=1,nbf
                  sum=zero
                  do 23 i=1,n
                     sum=sum+abs(grad(ioff+i,mu,coord))
   23             continue
                  gradbar(mu,coord)=sum/float(n)
                  gmax(coord)=max(gmax(coord),gradbar(mu,coord))
   24          continue
   25       continue
         endif
c
c        --- loop over spin components. first alpha, then beta.
         do 95 exch=1,2
c
c           --- multiply operator by weights into scratch area.
c               part common to all functionals.
            if(exch.eq.1) then
               call vmul(scr,fout(1+ioff,2),tmpgwt(1+ioff),n)
            else
               call vmul(scr,fout(1+ioff,4),tmpgwt(1+ioff),n)
            endif
c
c           --- get mean value of the operator on this block.
            sumo=zero
            do 30 i=1,n
               sumo=sumo+abs(scr(i))
   30       continue
            opbar=sumo/float(n)
c
c           --- loop over basis functions and generate XC matrix. 
            do 40 mu=1,nbf
               test=float(n)*opbar*phibar(mu)
               if(test*pmax.ge.kmcut) then
c                 --- form operator on phi(mu).
                  call vmul(tea(1,mu),scr,phi(1+ioff,mu),n)
                  do 35 nu=1,mu
                     kl=(mu*(mu-1)/2)+nu
c                    --- test the contribution of this block.
                     testnu=test*phibar(nu)
                     if (testnu.ge.kmcut) then
                        kmat(kl,exch)=kmat(kl,exch)
     $                               +sdot(n,tea(1,mu),1,
     $                                     phi(1+ioff,nu),1)
                     endif
   35             continue
               endif
   40       continue
c
c           --- perhaps do the gradient
            if(becke.or.lyp) then
c              --- form the operator
               if(exch.eq.1) then
                  call vadd(scr,fout(1+ioff,3),fout(1+ioff,3),n)
               else
                  call vadd(scr,fout(1+ioff,5),fout(1+ioff,5),n)
               endif
               do 90 coord=1,3
c                 --- form this component of the operator.
                  if(exch.eq.1) then
                     call vmul(queue(1,1),scr,dengrada(1+ioff,coord),n)
                     if(lyp) then
                        call vmul(scr2,fout(1+ioff,6),
     $                            dengradb(1+ioff,coord),n)
                        call vadd(queue(1,1),queue(1,1),scr2,n)
                     endif
                  else
                     call vmul(queue(1,1),scr,dengradb(1+ioff,coord),n)
                     if(lyp) then
                        call vmul(scr2,fout(1+ioff,6),
     $                            dengrada(1+ioff,coord),n)
                        call vadd(queue(1,1),queue(1,1),scr2,n)
                     endif
                  endif
c                 --- multiply by weights
                  call vmul(queue(1,1),queue(1,1),tmpgwt(1+ioff),n)
c                 --- get mean value of the operator on this block.
                  sumo=zero
                  do 50 i=1,n
                     sumo=sumo+abs(queue(i,1))
   50             continue
                  opbar=sumo/float(n)
c
c                 --- loop over basis functions and generate XC matrix. 
                  do 80 mu=1,nbf
c                    --- see if this block is negligible
                     test2=float(n)*opbar*phibar(mu)
                     test3=float(n)*opbar*gradbar(mu,coord)
                     test=test2*gmax(coord)+test3*pmax
                     if(test.ge.kmcut) then
c                       --- form operator on phi(mu) and gradphi(mu).
                        call vmul(queue(1,2),queue(1,1),
     $                            phi(1+ioff,mu),n)
                        call vmul(queue(1,3),queue(1,1),
     $                            grad(1+ioff,mu,coord),n)
                        do 70 nu=1,mu
                           kl=(mu*(mu-1)/2)+nu
c                          --- test the contribution of this block.
                           test4=test2*gradbar(nu,coord)
                           test5=test3*phibar(nu)
                           if (test4.ge.kmcut) then
                              kmat(kl,exch)=kmat(kl,exch)
     $                                     +sdot(n,queue(1,2),1,
     $                                        grad(1+ioff,nu,coord),1)
                           endif
                           if(test5.ge.kmcut) then
                              kmat(kl,exch)=kmat(kl,exch)
     $                                     +sdot(n,queue(1,3),1,
     $                                           phi(1+ioff,nu),1)
                           endif
   70                   continue
                     endif
   80             continue
   90          continue
            endif
   95    continue
c        --- increment the offset and do the next block.
         ioff=ioff+n
  100 continue
c
c
      return
      end
