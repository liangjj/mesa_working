*deck @(#)bfints.f	1.1 9/8/91
c***begin prologue     bfints
c***date written       xxxxxx   (yymmdd)
c***revision date      890427   (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            bound-free matrix elements
c***description        calculation of bound-free integrals
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       bfints
      subroutine bfints(hpvb,ovbf,grid,basis,hp,hd,ylm,vpot,cmat,rmat,
     1                  bmat,scrc,scrr,nlm,lch,mch,ngauss,ngch,npt,
     2                  nchan,nkept,list,lmax,nstri,dimlm,dimc,
     3                  dimbf,maxlm)
      implicit integer (a-z)
      real *8 hd, ylm, vpot, grid, basis, scrr, rmat
      complex *16  cmat, bmat, hpvb, scrc
      complex *16 ovbf, hp
      dimension ylm(npt,0:lmax,0:2*lmax), vpot(npt,nstri)
      dimension nlm(dimc), lch(dimlm,dimc), mch(dimlm,dimc)
      dimension cmat(maxlm,npt), bmat(npt,nkept), grid(4,npt)
      dimension rmat(maxlm,npt)
      dimension hpvb(1:maxlm,nkept,nchan,nchan), scrc(*)
      dimension ovbf(1:maxlm,nkept,nchan), scrr(*)
      dimension hp(npt,maxlm,nchan), hd(npt,maxlm,nchan)
      dimension basis(npt,nkept), list(dimbf)
      dimension ngauss(dimc), ngch(dimbf,dimc)
c----------------------------------------------------------------------c
c          first two loops over channel 1 and associated l, m          c
c
      do 10 ch1=1,nchan
         call czero(cmat,maxlm*npt)
         call filcf1(cmat,ylm,hp(1,1,ch1),nlm(ch1),lch(1,ch1),
     1               mch(1,ch1),maxlm,npt,lmax,dimlm)
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c             next two loops over channel 2 and associated             c
c                        bound orbitals                                c
         do 20 ch2=1,nchan
            call czero(bmat,npt*nkept)
            chlrg=max(ch1,ch2)
            chsml=ch1+ch2-chlrg
            itri=chlrg*(chlrg-1)/2+chsml
            do 30 bfn=1,ngauss(ch2)
               bfnch2 = ngch(bfn,ch2)
               blocte=list(bfnch2)
               if (blocte.eq.0) then
                   call lnkerr('bound orbital out of range')
               endif
c----------------------------------------------------------------------c
               do 40 grpt=1,npt
                  bmat(grpt,blocte)=vpot(grpt,itri) * basis(grpt,blocte)
   40          continue
   30       continue
            call cebc(scrc,cmat,bmat,nlm(ch1),npt,nkept)
            call cvadd(hpvb(1,1,ch1,ch2),maxlm,scrc,nlm(ch1),
     1                 nlm(ch1),nkept)
   20    continue
   10 continue
c----------------------------------------------------------------------c
c            compute matrix of (e - t)                                 c
c----------------------------------------------------------------------c
      do 70 ch1=1,nchan
         call czero(cmat,maxlm*npt)
         call rzero(rmat,maxlm*npt)
         call filc12(cmat,rmat,ylm,grid,hp(1,1,ch1),hd(1,1,ch1),
     1               nlm(ch1),lch(1,ch1),mch(1,ch1),maxlm,npt,lmax,
     2               dimlm)
         call czero(bmat,npt*nkept)
         do 80 bfn=1,ngauss(ch1)
            bfnch1=ngch(bfn,ch1)
            blocte=list(bfnch1)
            if (blocte.eq.0) then
                call lnkerr('bound orbital out of range')
            endif
            do 90 grpt=1,npt
               bmat(grpt,blocte)=basis(grpt,blocte)
   90       continue
   80    continue
         call cebc(scrc,cmat,bmat,nlm(ch1),npt,nkept)
         call cvadd(ovbf(1,1,ch1),maxlm,scrc,nlm(ch1),nlm(ch1),nkept)
         call ebc(scrr,rmat,basis,nlm(ch1),npt,nkept)
         call crvadd(hpvb(1,1,ch1,ch1),maxlm,scrr,nlm(ch1),nlm(ch1),
     1              nkept)
   70 continue
      return
      end








