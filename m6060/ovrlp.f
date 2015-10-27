*deck @(#)ovrlp.f
c***begin prologue     ovrlp
c***date written       xxxxxx   (yymmdd)
c***revision date      920409   (yymmdd)
c***keywords           kohn integrals
c***author             schneider, barry (nsf)
c***source             m6060
c***purpose            overlaps of bound and free functions
c***                             with vlamdas
c***references
c
c***routines called    iosys, util and mdutil
c***end prologue       ovrlp
      subroutine ovrlp(ovplm,ovmlm,ovlmb,ovpb,ovmb,basis,hp,ylm,
     1                 vlamda,grid,cmat,bmat,scrc,nlm,lch,mch,
     2                 ngauss,ngch,npt,nchan,nkept,list,lmax,dimlm,
     3                 dimc,dimbf,maxlm,nolam,bcond)
      implicit integer (a-z)
      character *(*) bcond
      real *8 ylm, basis, grid
      complex *16  cmat, bmat, scrc, vlamda, ovplm, ovmlm
      complex *16 ovlmb, hp, ovpb, ovmb
      dimension ylm(npt,0:lmax,0:2*lmax), vlamda(npt,nolam)
      dimension nlm(dimc), lch(dimlm,dimc), mch(dimlm,dimc)
      dimension cmat(maxlm,npt,2), bmat(npt,nkept), scrc(*)
      dimension ovplm(1:maxlm,nolam,nchan), ovmlm(1:maxlm,nolam,nchan)
      dimension ovpb(1:maxlm,nkept,nchan), ovmb(1:maxlm,nkept,nchan)
      dimension ovlmb(nolam,nkept,nchan), hp(npt,maxlm,nchan)
      dimension basis(npt,nkept), list(dimbf), ngauss(dimc)
      dimension ngch(dimbf,dimc), grid(4,npt)
      common /io/ inp, iout
c---------------------------------------------------------------------c
c             do overlap of free function with vlamda                 c
c---------------------------------------------------------------------c
      do 10 ch=1,nchan
         call czero(cmat,2*maxlm*npt)
         call filcf1(cmat(1,1,1),ylm,hp(1,1,ch),nlm(ch),lch(1,ch),
     1               mch(1,ch),maxlm,npt,lmax,dimlm)
         call fildf1(cmat(1,1,2),ylm,hp(1,1,ch),nlm(ch),lch(1,ch),
     1               mch(1,ch),maxlm,npt,lmax,dimlm,bcond)
         call cebc(scrc,cmat(1,1,1),vlamda,nlm(ch),npt,nolam)
         call cvadd(ovplm(1,1,ch),maxlm,scrc,nlm(ch),nlm(ch),nolam)
         call cebc(scrc,cmat(1,1,2),vlamda,nlm(ch),npt,nolam)
         call cvadd(ovmlm(1,1,ch),maxlm,scrc,nlm(ch),nlm(ch),nolam)
   10 continue
c---------------------------------------------------------------------c
c             do overlap of bound function with vlamda                c
c---------------------------------------------------------------------c
      do 20 ch=1,nchan
         call czero(bmat,npt*nkept)
         do 30 bfn=1,ngauss(ch)
            bfnch=ngch(bfn,ch)
            blocte=list(bfnch)
            if(blocte.eq.0) then
               call lnkerr('bound orbital out of range')
            endif
            do 40 grpt=1,npt
               bmat(grpt,blocte)=basis(grpt,blocte)
   40       continue
   30    continue
         call cebtc(scrc,vlamda,bmat,nolam,npt,nkept)
         call cvadd(ovlmb(1,1,ch),nolam,scrc,nolam,nolam,nkept)
   20 continue
c---------------------------------------------------------------------c
c               do overlap of bound and free functions                c
c---------------------------------------------------------------------c
      do 50 ch=1,nchan
         call czero(cmat,2*maxlm*npt)
         call filcf1(cmat(1,1,1),ylm,hp(1,1,ch),nlm(ch),lch(1,ch),
     1               mch(1,ch),maxlm,npt,lmax,dimlm)
         call fildf1(cmat(1,1,2),ylm,hp(1,1,ch),nlm(ch),lch(1,ch),
     1               mch(1,ch),maxlm,npt,lmax,dimlm,bcond)
         call czero(bmat,npt*nkept)
         do 60 bfn=1,ngauss(ch)
            bfnch=ngch(bfn,ch)
            blocte=list(bfnch)
            if(blocte.eq.0) then
               call lnkerr('bound orbital out of range')
            endif
            do 70 grpt=1,npt
               bmat(grpt,blocte)=grid(4,grpt)*basis(grpt,blocte)
   70       continue
   60    continue
         call cebc(scrc,cmat(1,1,1),bmat,nlm(ch),npt,nkept)
         call cvadd(ovpb(1,1,ch),maxlm,scrc,nlm(ch),nlm(ch),nkept)
         call cebc(scrc,cmat(1,1,2),bmat,nlm(ch),npt,nkept)
         call cvadd(ovmb(1,1,ch),maxlm,scrc,nlm(ch),nlm(ch),nkept)
   50 continue
      return
      end





