*deck @(#)cint.f	5.1  11/6/94
      subroutine cint(t,ti,b,bi,tgrad,wt,iptci,nmix,lmixm,
     1                nspinm)
      implicit real*8(a-h,o-z)
c
      integer iptci(2,10,*)
      integer lmixm(*)
      real*8 bi(*),wt(10,*),t(nmix),b(nmix)
      real*8 tgrad(nmix,*),ti(*)
      real*8 sdot
c
      ix=0
      kx=0
      do 11 k=1,nspinm
         lmix=lmixm(k)
         if(lmix.ne.0) then
            do 10 i=1,lmix
c
               ic=iptci(1,i,k)
               jc=iptci(2,i,k)
c
               kx=kx+1
c
c              phase changed nov. 1987 bhl
c
               xfac=-(wt(jc,k)-wt(ic,k))*bi(i+ix)
c
               call saxpy(nmix,xfac,tgrad(1,kx),1,t,1)
c
               xr=sdot(nmix,b,1,tgrad(1,kx),1)
c
c              phase changed nov. 1987 bhl
c
               ti(i+ix)=ti(i+ix)-(wt(jc,k)-wt(ic,k))*xr
  10        continue
         endif
c
         ix=ix+lmix
  11  continue
c
      return
      end
