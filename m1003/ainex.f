*deck @(#)ainex.f	5.1  11/6/94
      subroutine ainex(b,bi,t,ti,thc,wt,iptci,ncsfm,nwtm,lmixm,
     $                 nspinm)
      implicit real*8(a-h,o-z)
c
      integer iptci(2,10,*),lmixm(*)
      integer ncsfm(*),nwtm(*)
      real*8 bi(*),ti(*),wt(10,*)
      real*8 b(*),t(*),thc(*)
      real*8 sdot
c
      ix=0
      ib=1
c
      do 11 k=1,nspinm
         ncsf=ncsfm(k)
         nwt=nwtm(k)
         lmix=lmixm(k)
         if(lmix.ne.0) then
            do 10 i=1,lmix
c
               ic=iptci(1,i,k)
               jc=iptci(2,i,k)
               ipt=(ic-1)*ncsf+ib
               jpt=(jc-1)*ncsf+ib
c
               rj=sdot(ncsf,b(jpt),1,thc(ipt),1)
               ri=sdot(ncsf,b(ipt),1,thc(jpt),1)
c
               ti(ix+i)=ti(ix+i)+wt(ic,k)*ri-wt(jc,k)*rj
c
               xi=wt(ic,k)*bi(i+ix)
               xj=-wt(jc,k)*bi(i+ix)
c
               call saxpy(ncsf,xi,thc(jpt),1,t(ipt),1)
               call saxpy(ncsf,xj,thc(ipt),1,t(jpt),1)
c
  10        continue
         endif
c
         ib=ib+ncsf*nwt
         ix=ix+lmix
c
  11  continue
c
      return
      end
