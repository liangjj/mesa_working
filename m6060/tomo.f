*deck @(#)tomo.f	1.1 9/8/91
c***begin prologue     tomo
c***date written       890529   (yymmdd)
c***revision date               (yymmdd)
c***keywords           transformation, kohn
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            bound-free ao to mo transformation
c***
c***references         none
c
c***routines called    cebc
c***end prologue       tomo
 
      subroutine tomo(trans,ovlmb,ovlmbm,ovpb,ovpbm,ovmb,ovmbm,
     1                scrc,scrcc,maxlm,nolam,ncon,nmo,nchan,ngauss,
     2                ngch,nkept,list,dimbf,dimc)
      implicit integer (a-z)
      real *8 trans
      complex *16 ovlmb, ovlmbm, ovpb, ovmb, ovpbm, ovmbm, scrc
      complex *16 scrcc
      dimension ovlmb(nolam,nkept,nchan), ovlmbm(nolam,nmo,nchan)
      dimension ovpb(1:maxlm,nkept,nchan), ovmb(1:maxlm,nkept,nchan)
      dimension ovpbm(1:maxlm,nmo,nchan), ovmbm(1:maxlm,nmo,nchan) 
      dimension trans(ncon,nmo), list(dimbf), scrc(nolam,ncon)
      dimension ngauss(dimc), ngch(dimbf,dimc)
      dimension scrcc(1:maxlm,ncon)
      do 10 ch=1,nchan
         call czero(scrc,nolam*ncon)
         do 20 bfn=1,ngauss(ch)
            bfnch = ngch(bfn,ch)
            blocte=list(bfnch)
            do 30 lam=1,nolam
               scrc(lam,bfnch)=ovlmb(lam,blocte,ch)
   30       continue
   20    continue
         call ecbcx(ovlmbm(1,1,ch),nolam,scrc,nolam,trans,
     1              ncon,nolam,ncon,nmo)
   10 continue
      do 50 ch=1,nchan
         call czero(scrcc,nolam*ncon)
         do 60 bfn=1,ngauss(ch)
            bfnch = ngch(bfn,ch)
            blocte=list(bfnch)
            do 70 lm=1,maxlm
               scrcc(lm,bfnch)=ovpb(lm,blocte,ch)
   70       continue
   60    continue
         call ecbcx(ovpbm(1,1,ch),maxlm,scrcc,maxlm,trans,
     1              ncon,maxlm,ncon,nmo)
   50 continue
      do 90 ch=1,nchan
         call czero(scrcc,maxlm*ncon)
         do 100 bfn=1,ngauss(ch)
            bfnch = ngch(bfn,ch)
            blocte=list(bfnch)
            do 110 lm=1,maxlm
               scrcc(lm,bfnch)=ovmb(lm,blocte,ch)
  110       continue
  100    continue
         call ecbcx(ovmbm(1,1,ch),maxlm,scrcc,maxlm,trans,
     1              ncon,maxlm,ncon,nmo)
   90 continue
      return
      end



