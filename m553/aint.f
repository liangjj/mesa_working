*deck @(#)aint.f	5.1  11/6/94
      subroutine aint(ti,bi,wt,wener,iptci,lmixm,nwtm,nspinm)
      implicit real*8(a-h,o-z)
      real*8 ti(*),bi(*),wt(10,*),wener(10,*)
      real*8 small
      integer iptci(2,10,*),nwtm(*)
      integer lmixm(*)
      parameter (small=1.d-3)
c
c fix for degenerate weights
c
      ix=0
      do 11 k=1,nspinm
         lmix=lmixm(k)
         if(lmix.ne.0) then
            do 10 i=1,lmix
               ic=iptci(1,i,k)
               jc=iptci(2,i,k)
               if(abs(wt(ic,k)-wt(jc,k)).gt.small) then
                  ti(i+ix)=ti(i+ix)+
     $            bi(i+ix)*(wt(ic,k)-wt(jc,k))*(wener(jc,k)-wener(ic,k))
               else
                  ti(i+ix)=ti(i+ix)+bi(i+ix)
               end if
  10        continue
         endif
         ix=ix+lmix
  11  continue
c
      return
      end
