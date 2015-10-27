*deck @(#)ffint1.f	1.1 9/8/91
c***begin prologue     ffint1
c***date written       890529   (yymmdd)
c***revision date               (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            free-free matrix elements of kinetic energy
c***description        calculation of free-free integrals
c***references
c***
c***
c***routines called    iosys, util and mdutil
c***end prologue       ffint1
      subroutine ffint1(ovpp,ovpm,hpvhp,grid,hp,hd,ylm,cmatt,rmat,
     1                  rtrick,dmat,nlm,lch,mch,npt,nchan,lmax,dimlm,
     2                  dimc,maxlm)
      implicit integer (a-z)
      common /io/ inp,iout 
      real *8 ylm, grid, hd, sdot, rmat, rtrick, dmat
      complex *16 cmatt, hpvhp, ovpp, ovpm
      complex *16 hp, ai, cdotu
      dimension ylm(npt,0:lmax,0:2*lmax), hd(npt,maxlm,nchan)
      dimension nlm(dimc), lch(dimlm,dimc), mch(dimlm,dimc)
      dimension hpvhp(1:maxlm,nchan), ovpp(1:maxlm,nchan)
      dimension cmatt(npt,maxlm), grid(4,npt), dmat(npt,maxlm)
      dimension ovpm(1:maxlm,nchan), hp(npt,maxlm,nchan)
      dimension rmat(maxlm,2), rtrick(2*npt,maxlm)
      ai=dcmplx(0.d+00,1.d+00)
c----------------------------------------------------------------------c
c            compute matrix of (kchan**2/2 - t)                        c
c            diagonal in channel indices so its simple                 c
c       rtrick (real) and cmatt (complex) must be equivalenced         c
c       in the call to this routine for it to function properly        c
c----------------------------------------------------------------------c
      do 10 ch1=1,nchan
         clst=ch1*(ch1+1)/2
         call czero(cmatt,maxlm*npt)
         do 20 nolm1=1,nlm(ch1)
            l1=lch(nolm1,ch1)
            m1=mch(nolm1,ch1)
            do 30 grpt=1,npt
               cmatt(grpt,nolm1) = grid(4,grpt)*ylm(grpt,l1,m1)
     1                             *ylm(grpt,l1,m1)
     2                             *hp(grpt,nolm1,ch1)
               dmat(grpt,nolm1)=imag(hp(grpt,nolm1,ch1))
   30       continue
   20    continue
         do 40 nolm1=1,nlm(ch1)
            rmat(nolm1,1)=sdot(npt,rtrick(1,nolm1),2,hd(1,nolm1,ch1),1)
            rmat(nolm1,2)=sdot(npt,rtrick(2,nolm1),2,hd(1,nolm1,ch1),1)
            hpvhp(nolm1,ch1)=hpvhp(nolm1,ch1) + rmat(nolm1,1) 
     1                                        + ai*rmat(nolm1,2)
            ovpp(nolm1,ch1)=ovpp(nolm1,ch1)+cdotu(npt,cmatt(1,nolm1),1,
     1                                            hp(1,nolm1,ch1),1)
   40    continue
         do 50 nolm1=1,nlm(ch1)
            rmat(nolm1,1)=sdot(npt,rtrick(1,nolm1),2,dmat(1,nolm1),1)
            rmat(nolm1,2)=sdot(npt,rtrick(2,nolm1),2,dmat(1,nolm1),1)
            ovpm(nolm1,ch1)=ovpm(nolm1,ch1) + rmat(nolm1,1)
     1                                      +ai*rmat(nolm1,2)
   50    continue
   10 continue
      return
      end








