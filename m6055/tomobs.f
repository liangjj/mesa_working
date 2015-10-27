*deck @(#)tomobs.f	1.1 9/8/91
c***begin prologue     tomobs
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
c***end prologue       tomobs
 
      subroutine tomobs(trans,ovbf,ovbm,hpvb,hpvm,scrc,maxlm,ncon,nmo,
     1                  nchan,ngauss,ngch,nkept,list,dimbf,dimc,exdim)
      implicit integer (a-z)
      real *8 trans
      complex *16 ovbf, ovbm, hpvb, hpvm, scrc
      dimension ovbf(1:maxlm,nkept,nchan), ovbm(1:maxlm,nmo,nchan)
      dimension hpvb(1:maxlm,nkept,nchan,exdim), trans(ncon,nmo)
      dimension hpvm(1:maxlm,nmo,nchan,exdim), list(dimbf)
      dimension scrc(1:maxlm,ncon,nchan,exdim)
      dimension ngauss(dimc), ngch(dimbf,dimc)
      nzero=maxlm*ncon*nchan*exdim
      call czero(scrc,nzero)
      do 10 ch1=1,nchan
         do 20 bfn=1,ngauss(ch1)
            bfnch = ngch(bfn,ch1)
            blocte=list(bfnch)
            do 30 lm=1,maxlm
               scrc(lm,bfnch,ch1,1)=ovbf(lm,blocte,ch1)
   30       continue
   20    continue
   10 continue
      do 40 ch1=1,nchan
         call ecbcx(ovbm(1,1,ch1),maxlm,scrc(1,1,ch1,1),maxlm,trans,
     1             ncon,maxlm,ncon,nmo)
   40 continue
      call czero(scrc,nzero)
      do 50 ch1=1,nchan
         do 60 ch2=1,exdim
            do 70 bfn=1,ngauss(ch2)
               bfnch = ngch(bfn,ch2)
               blocte=list(bfnch)
               do 80 lm=1,maxlm
                  scrc(lm,bfnch,ch1,ch2)=hpvb(lm,blocte,ch1,ch2)
   80          continue
   70       continue
   60    continue
   50 continue
      do 90 ch1=1,nchan
         do 100 ch2=1,exdim
            call ecbcx(hpvm(1,1,ch1,ch2),maxlm,scrc(1,1,ch1,ch2),maxlm,
     1                 trans,ncon,maxlm,ncon,nmo)
  100    continue
   90 continue
      return
      end
