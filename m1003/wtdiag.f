*deck @(#)wtdiag.f	5.1  11/6/94
      subroutine wtdiag(b,diag,di,thc,c0,wt,wener,iptci,
     1                  nwtm,ncsfm,lmixm,nspinm)
      implicit real*8(a-h,o-z)
c
      dimension b(*),diag(*),di(*),wt(10,*),iptci(2,10,*)
      dimension wener(10,*),thc(*),c0(*),ncsfm(*),nwtm(*),lmixm(*)
      real*8 sdot
c
      small=1.d-08
c
      ix=1
      jx=1
      mx=0
      do 101 k=1,nspinm
         ncsf=ncsfm(k)
         nwt=nwtm(k)
         do 100 i=1,nwt
            call vmove(diag(ix),b(jx),ncsf)
            call sadd(diag(ix),diag(ix),-wener(i,k),ncsf)
            ix=ix+ncsf
  100    continue
         jx=jx+ncsf
  101 continue

c
c
      kx=1
      do 51 k=1,nspinm
         ncsf=ncsfm(k)
         nwt=nwtm(k)
         jx=kx
         do 50 i=1,nwt
            xfac=2.d0*wener(i,k)
            if(wt(i,k).lt.small) go to 21
            ix=kx
            do 10 n=1,nwt
               if(wt(n,k).lt.small) go to 10
               rr=sdot(ncsf,thc(ix),1,c0(ix),1)
c              write(iout,10005) wener(i,k),rr
c10005         format(' wener rr ',2(2x,f12.8))
c              call prntv(c0(ix),ncsf)
c              call prntv(thc(ix),ncsf)
               call vsamul(diag(jx),thc(ix),c0(ix),-2.d0,ncsf)
               call vsamul(diag(jx),c0(ix),c0(ix),xfac,ncsf)
  10        ix=ix+ncsf
c
            ix=kx
            do 20 n=1,nwt
               if(i.eq.n) go to 88009
               if(wt(n,k).lt.small) go to 88009
               xfac=wener(n,k)-wener(i,k)
               call vsamul(diag(jx),c0(ix),c0(ix),xfac,ncsf)
88009          ix=ix+ncsf
  20        continue
  21        continue
c
            if(wt(i,k).lt.small) go to 41
            call sscal(ncsf,wt(i,k),diag(jx),1)
            go to 50
  41        continue
            do 42 j=1,ncsf
               diag(jx-1+j)=1.d0
  42        continue
c
  50     jx=jx+ncsf
cc
cc
         lmix=lmixm(k)
         if(lmix.ne.0) then
            do 60 m=1,lmix
               ic=iptci(1,m,k)
               jc=iptci(2,m,k)
               di(m+mx)=(wt(ic,k)-wt(jc,k))*(wener(jc,k)-wener(ic,k))
  60        continue
         endif
         kx=kx+ncsf*nwt
         mx=mx+lmix
  51  continue
c
      return
      end
