*deck @(#)bfint1.f	1.1 9/8/91
c***begin prologue     bfint1
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            bound-free matrix elements of kinetic energy
c***description        calculation of bound-free integrals
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       bfint1
      subroutine bfint1(hpvb,ovbf,grid,basis,hp,hd,ylm,cmat,rmat,
     1                  bmat,scrc,scrr,nlm,lch,mch,ngauss,ngch,npt,
     2                  nchan,nkept,list,lmax,dimlm,dimc,dimbf,maxlm)
      implicit integer (a-z)
      real *8 hd, ylm, grid, basis, scrr, rmat
      complex *16  cmat, bmat, hpvb, scrc
      complex *16 ovbf, hp
      dimension ylm(npt,0:lmax,0:2*lmax)
      dimension nlm(dimc), lch(dimlm,dimc), mch(dimlm,dimc)
      dimension cmat(maxlm,npt), bmat(npt,nkept), grid(4,npt)
      dimension rmat(maxlm,npt)
      dimension hpvb(1:maxlm,nkept,nchan), scrc(*)
      dimension ovbf(1:maxlm,nkept,nchan), scrr(*)
      dimension hp(npt,maxlm,nchan), hd(npt,maxlm,nchan)
      dimension basis(npt,nkept), list(dimbf)
      dimension ngauss(dimc), ngch(dimbf,dimc)
c----------------------------------------------------------------------c
c            compute matrix of (kchan**2/2 - t)                        c
c----------------------------------------------------------------------c
      do 10 ch1=1,nchan
         call czero(cmat,maxlm*npt)
         call rzero(rmat,maxlm*npt)
         call filc12(cmat,rmat,ylm,grid,hp(1,1,ch1),hd(1,1,ch1),
     1               nlm(ch1),lch(1,ch1),mch(1,ch1),maxlm,npt,lmax,
     2               dimlm)
         call czero(bmat,npt*nkept)
         do 20 bfn=1,ngauss(ch1)
            bfnch1=ngch(bfn,ch1)
            blocte=list(bfnch1)
            if (blocte.eq.0) then
                call lnkerr('bound orbital out of range')
            endif
            do 30 grpt=1,npt
               bmat(grpt,blocte)=basis(grpt,blocte)
   30       continue
   20    continue
         call cebc(scrc,cmat,bmat,nlm(ch1),npt,nkept)
         call cvadd(ovbf(1,1,ch1),maxlm,scrc,nlm(ch1),nlm(ch1),nkept)
         call ebc(scrr,rmat,basis,nlm(ch1),npt,nkept)
         call crvadd(hpvb(1,1,ch1),maxlm,scrr,nlm(ch1),nlm(ch1),
     1               nkept)
   10 continue
      return
      end








