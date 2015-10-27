*deck errfac.f
c***begin prologue     errfac
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             tstdiis
c***purpose            calculate the scalar products needed to form
c***                   the error matrix on the fly.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       errfac
      subroutine errfac(psi,hpsi,hh,pp,ph,hp,n,iter,maxit)
      implicit integer (a-z)
      real*8 psi, hpsi, hh, pp, ph, hp, sdot
      dimension psi(n,*), hpsi(n,*), hh(maxit,maxit)
      dimension pp(maxit,maxit), ph(maxit,maxit), hp(maxit,maxit)
      common/io/inp, iout
      do 10 i=1,iter-1
         hh(i,iter) = sdot(n,hpsi(1,i),1,hpsi(1,iter),1)
         pp(i,iter) = sdot(n,psi(1,i),1,psi(1,iter),1)
         ph(i,iter) = sdot(n,psi(1,i),1,hpsi(1,iter),1)
         hp(i,iter) = sdot(n,hpsi(1,i),1,psi(1,iter),1)
         hh(iter,i) = hh(i,iter)
         pp(iter,i) = pp(i,iter)
         ph(iter,i) = ph(i,iter)
         hp(iter,i) = hp(i,iter)
 10   continue
         hh(iter,iter) = sdot(n,hpsi(1,iter),1,hpsi(1,iter),1)
         pp(iter,iter) = sdot(n,psi(1,iter),1,psi(1,iter),1)
         ph(iter,iter) = sdot(n,psi(1,iter),1,hpsi(1,iter),1)
         hp(iter,iter) = sdot(n,hpsi(1,iter),1,psi(1,iter),1)         
      return 
      end       
