*deck @(#)ffints.f	1.1 9/8/91
c***begin prologue     ffints
c***date written       890529   (yymmdd)
c***revision date               (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            free-free matrix elements
c***description        calculation of free-free integrals
c***references
c***
c***
c***routines called    iosys, util and mdutil
c***end prologue       ffints
      subroutine ffints(hpvhp,hpvhm,grid,hp,hd,ylm,vpot,cmat,cmatt,
     1                  dmat,rmat,rtrick,scrc,nlm,lch,mch,npt,nchan,
     2                  lmax,nstri,dimlm,dimc,maxlm,bcond)
      implicit integer (a-z)
      common /io/ inp,iout 
      character *(*) bcond
      real *8 ylm, vpot, grid, hd, sdot, rmat, rtrick
      complex *16 cmat, cmatt, dmat, hpvhp, scrc
      complex *16 hpvhm, hp, ai
      dimension ylm(npt,0:lmax,0:2*lmax), vpot(npt,nstri)
      dimension nlm(dimc), lch(dimlm,dimc), mch(dimlm,dimc)
      dimension cmat(maxlm,npt,2), dmat(maxlm,npt,2), grid(4,npt)
      dimension hpvhp(1:maxlm,1:maxlm,nstri), scrc(maxlm,maxlm)
      dimension hpvhm(1:maxlm,1:maxlm,nchan,nchan), cmatt(npt,maxlm)
      dimension hp(npt,maxlm,nchan), hd(npt,maxlm,nchan)
      dimension rmat(maxlm,2), rtrick(2*npt,maxlm)
      ai=dcmplx(0.d+00,1.d+00)
c
c
c compute potential matrices
c
c----------------------------------------------------------------------c
c         first two loops over channel 1 and associated l, m           c
c         to load some temporary matrices for vectorization            c
      do 10 ch1=1,nchan
         call czero(cmat(1,1,1),maxlm*npt)
         call czero(dmat(1,1,1),maxlm*npt)
         call filcf1(cmat(1,1,1),ylm,hp(1,1,ch1),nlm(ch1),
     1               lch(1,ch1),mch(1,ch1),maxlm,npt,lmax,dimlm)
         call fildf1(dmat(1,1,1),ylm,hp(1,1,ch1),nlm(ch1),
     1               lch(1,ch1),mch(1,ch1),maxlm,npt,lmax,dimlm,
     2               bcond)
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c          next two loops over channel 2 and associated l, m           c
c          same strategy on temporary matrices                         c
         do 20 ch2=1,ch1
            call czero(cmat(1,1,2),maxlm*npt)
            call czero(dmat(1,1,2),maxlm*npt)
            ist = ch1*(ch1-1)/2 + ch2
            call filcf2(cmat(1,1,2),ylm,hp(1,1,ch2),vpot(1,ist),
     1	                nlm(ch2),lch(1,ch2),mch(1,ch2),maxlm,npt,
     2                  lmax,dimlm)
            call fildf2(dmat(1,1,2),ylm,hp(1,1,ch2),vpot(1,ist),
     1                  nlm(ch2),lch(1,ch2),mch(1,ch2),maxlm,npt,
     2                  lmax,dimlm,bcond)
c----------------------------------------------------------------------c
c                  accumulate into matrix elements                     c
c----------------------------------------------------------------------c
            call cebct(scrc,cmat(1,1,1),cmat(1,1,2),nlm(ch1),
     1                 npt,nlm(ch2))
            call cvadd(hpvhp(1,1,ist),maxlm,scrc,nlm(ch1),nlm(ch1),
     1                 nlm(ch2))
            if (ch1.ne.ch2) then
                call cebct(scrc,cmat(1,1,1),dmat(1,1,2),nlm(ch1),
     1                     npt,nlm(ch2))
                call cvadd(hpvhm(1,1,ch1,ch2),maxlm,scrc,nlm(ch1),
     1                     nlm(ch1),nlm(ch2))
            endif
            call cebct(scrc,cmat(1,1,2),dmat(1,1,1),nlm(ch2),npt,
     1                 nlm(ch1))
            call cvadd(hpvhm(1,1,ch2,ch1),maxlm,scrc,nlm(ch2),
     1                 nlm(ch2),nlm(ch1))
   20    continue
   10 continue
c----------------------------------------------------------------------c
c            compute matrix of (kchan**2/2 - t)                        c
c            diagonal in channel indices so its simplier than above    c
c       rtrick (real) and cmatt (complex) must be equivalenced         c
c       in the call to this routine for it to function properly        c
c----------------------------------------------------------------------c
      do 30 ch1=1,nchan
         clst=ch1*(ch1+1)/2
         call czero(cmatt,maxlm*npt)
         do 40 nolm1=1,nlm(ch1)
            l1=lch(nolm1,ch1)
            m1=mch(nolm1,ch1)
            do 50 grpt=1,npt
               cmatt(grpt,nolm1) = grid(4,grpt)*ylm(grpt,l1,m1)
     1                             *ylm(grpt,l1,m1)
     2                             *hp(grpt,nolm1,ch1)
   50       continue
   40    continue
         do 60 nolm1=1,nlm(ch1)
            rmat(nolm1,1)=sdot(npt,rtrick(1,nolm1),2,
     1                                hd(1,nolm1,ch1),1)
            rmat(nolm1,2)=sdot(npt,rtrick(2,nolm1),2,
     1                             hd(1,nolm1,ch1),1)
            hpvhp(nolm1,nolm1,clst)=hpvhp(nolm1,nolm1,clst) +
     1                               rmat(nolm1,1) + ai*rmat(nolm1,2)
   60    continue
         if(bcond.eq.'s-matrix') then
            do 70 nolm1=1,nlm(ch1)
               hpvhm(nolm1,nolm1,ch1,ch1)=hpvhm(nolm1,nolm1,ch1,ch1) +
     1                                    rmat(nolm1,1) + 
     2                                              ai*rmat(nolm1,2)
   70       continue
         endif
   30 continue
      return
      end
