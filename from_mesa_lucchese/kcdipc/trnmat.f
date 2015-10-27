      subroutine trnmat(hpvbx,hpvby,hpvbz,hpvbtx,hpvbty,hpvbtz,
     $hpvb,hpvbt,hbx,hby,hbz,hbtx,tr,nchan2,nlm,
     $ngauss,nmotot,nchnl,nbfmax,lm2,iene)
      real*8 hpvbx(lm2,nbfmax,nchnl),hpvby(lm2,nbfmax,nchnl)
     $,hpvbz(lm2,nbfmax,nchnl),hpvbtx(lm2,nbfmax,nchnl)
     $,hpvbty(lm2,nbfmax,nchnl),hpvbtz(lm2,nbfmax,nchnl)
     $,hpvb(lm2,nbfmax,nchnl),hpvbt(lm2,nbfmax,nchnl)
      real*8 hbx(nbfmax,nbfmax),hby(nbfmax,nbfmax),
     $hbz(nbfmax,nbfmax),hbtx(nbfmax,nbfmax),tr(nbfmax,nbfmax)
      integer nlm(nchnl)
      do 1 ic=1,nchan2
      nl=nlm(ic)
      call mxma(hpvbx(1,1,ic),2,lm2,tr,1,nbfmax,hpvbtx(1,1,ic),
     $          2,lm2,nl,ngauss,nmotot)
      call mxma(hpvbx(2,1,ic),2,lm2,tr,1,nbfmax,hpvbtx(2,1,ic),
     $          2,lm2,nl,ngauss,nmotot)
      call mxma(hpvby(1,1,ic),2,lm2,tr,1,nbfmax,hpvbty(1,1,ic),
     $          2,lm2,nl,ngauss,nmotot)
      call mxma(hpvby(2,1,ic),2,lm2,tr,1,nbfmax,hpvbty(2,1,ic),
     $          2,lm2,nl,ngauss,nmotot)
      call mxma(hpvbz(1,1,ic),2,lm2,tr,1,nbfmax,hpvbtz(1,1,ic),
     $          2,lm2,nl,ngauss,nmotot)
      call mxma(hpvbz(2,1,ic),2,lm2,tr,1,nbfmax,hpvbtz(2,1,ic),
     $          2,lm2,nl,ngauss,nmotot)
      call mxma(hpvb(1,1,ic),2,lm2,tr,1,nbfmax,hpvbt(1,1,ic),
     $          2,lm2,nl,ngauss,nmotot)
      call mxma(hpvb(2,1,ic),2,lm2,tr,1,nbfmax,hpvbt(2,1,ic),
     $          2,lm2,nl,ngauss,nmotot)
 1    continue
      if(iene.eq.1)then
      call mxma(hbx,1,nbfmax,tr,1,nbfmax,hbtx,1,nbfmax,ngauss,
     $          ngauss,nmotot)
      call mxma(tr,nbfmax,1,hbtx,1,nbfmax,hbx,1,nbfmax,nmotot,
     $          ngauss,nmotot)
      call mxma(hby,1,nbfmax,tr,1,nbfmax,hbtx,1,nbfmax,ngauss,
     $          ngauss,nmotot)
      call mxma(tr,nbfmax,1,hbtx,1,nbfmax,hby,1,nbfmax,nmotot,
     $          ngauss,nmotot)
      call mxma(hbz,1,nbfmax,tr,1,nbfmax,hbtx,1,nbfmax,ngauss,
     $          ngauss,nmotot)
      call mxma(tr,nbfmax,1,hbtx,1,nbfmax,hbz,1,nbfmax,nmotot,
     $          ngauss,nmotot)
      endif
      return
      end
