c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{PREPIT: Prepare for Davidson Code}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck prepit.f
c***begin prologue     prepit
c***date written       010828   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           iterative eigenvalue
c***author             schneider, barry (nsf)
c***source             
c***purpose            prepare for eigdvr routine
c***references         
c
c***routines called    
c***end prologue       prepit
      subroutine prepit(pham,px,pv,dim,n,nd,nroot,ntrial,nattim,
     1                  cnverg,thresh,precon,maxit,maxvec,nblck,dvdprt)
      implicit integer (a-z)
      character*8 key
      real*8 hx, hy, hz, vpot, trial, h, dvd
      real*8 cnverg, thresh
#ifdef DECPOINTER
      integer*8 pham, pv, px
      integer*8 phx, phy, phz, pre
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 pham, pv, px
      integer*4 phx, phy, phz, pre
#endif SGIPOINTER
      logical dvdprt, zeroit, prnt
      character*(*) precon
      dimension nd(dim+1), pham(*), pv(dim+1), px(*) 
      dimension dvdprt(12)
      common/io/inp, iout
      pointer (ptrial,trial(1)), (ptrial,itrial(1))
      pointer (pdvd,dvd(1))
      pointer (phx,hx(1))
      pointer (phy,hy(1))
      pointer (phz,hz(1))
      pointer (pvpot,vpot(1))
      pointer (pre,h(1))
c
c     calculate the perturbation potential
c
      key='$vpert'
      zeroit=.true.
      call vpert(px,pv,nd,dim,key,vword,zeroit,prnt)
      pvpot=pv(dim+1)
c
      phx=pham(1)
      phy=pham(2)
      phz=pham(3)
      vec=1
      eigv=vec+n*ntrial
      need=wpadti(eigv+ntrial) 

      if(precon.ne.'block') then
c
c        prepare for preconditioning and trial vectors based on the separable
c        hamiltonian ( hx + hy + hz )
c
         call seppre(pre,hx,hy,hz,nd(1),nd(2),nd(3),n,dim)
         if(dim.ge.2) then
            ind=need
            need=need+n*dim
         endif
         call getmem(need,ptrial,ngot,'trial',0) 
         if(dim.eq.2) then
            call setnd2(itrial(ind),nd(1),nd(2),n)
         endif
         if(dim.eq.3) then
            call setnd3(itrial(ind),nd(1),nd(2),nd(3),n)
         endif         
         u1=1
         eig1=u1+nd(1)*nd(1)
         u2=eig1+nd(1)
         eig2=u2+nd(2)*nd(2)
         u3=eig2+nd(2)
         eig3=u3+nd(3)*nd(3)
         call vtrsep(h(u1),h(u2),h(u3),h(eig1),h(eig2),h(eig3),
     1               trial(vec),trial(eigv),itrial(ind),
     2               nd(1),nd(2),nd(3),dim,n,ntrial,.false.) 
c
      elseif(precon.eq.'block') then
c
c        prepare for preconditioning based on block diagonalization
c        of full hamiltonian.
c
         call blkpre(pre,hx,hy,hz,vpot,nd(1),nd(2),nd(3),n,nblck,dim)
         u=1
         eig=u+nblck*nblck
         ind=need
         need=need+n
         call getmem(need,ptrial,ngot,'trial',0) 
         call vtrblk(h(eig),h(u),trial(eigv),trial(vec),itrial(ind),
     1               n,ntrial,nblck,.false.)
      endif
c
c     the trial vectors are stored, we can get rid of the memory.
c
      call getmem(-ngot,ptrial,idum,'trial',idum)
c
c     now set up the memory for the davidson subroutine.
c
      diag=1
      eig=diag+n
      vec=eig+nroot
      hvec=vec+n*maxvec
      resid=hvec+n*maxvec
      b=resid+n*maxvec
      bwrk=b+maxvec*maxvec
      eigwrk=bwrk+maxvec*maxvec
      work=eigwrk+maxvec
      svec=work+max(5*maxvec,n*maxvec)
      need=wpadti(svec+maxvec*maxvec)
      call getmem(need,pdvd,nword,'dvd',0)
c
c     store the full diagonal element and then zero
c     the diagonal elements of the one-dimensional matrices.
c
      call diagnl(dvd(diag),hx,hy,hz,vpot,n,nd(1),nd(2),nd(3),dim)
c
c     the potential is no longer needed, its on the diagonal.
c
      call getmem(-vword,pvpot,ndum,'vint',ndum)
c
      call dvrdvd(hx,hy,hz,dvd(diag),dvd(eig),dvd(vec),dvd(hvec),
     1            dvd(resid),dvd(b),dvd(bwrk),dvd(eigwrk),dvd(work),
     2            dvd(svec),cnverg,thresh,dim,nd(1),nd(2),nd(3),
     3            n,nroot,ntrial,nattim,maxit,maxvec,nblck,
     4            dvdprt,precon,pre)
      call getmem(-nword,pdvd,idum,'dvd',idum)
      return
      end       


