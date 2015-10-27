*deck @(#)vcalcr.f	5.1  11/6/94
      function vcalcr(nsite,zm,zl,site,vo)
      implicit none
c
      real*8 vcalcr
      integer nsite
      real*8 zm(nsite),zl(3,nsite),site(3,nsite),vo(4)
      integer isite
      real*8 x,y,z,dsqr,d1,d3
      real*8 zero,one
c
      parameter (zero=0.0d+00,one=1.0d+00)
c
c
      vcalcr=zero
      do 10 isite=1,nsite
         x=vo(1)-site(1,isite)
         y=vo(2)-site(2,isite)
         z=vo(3)-site(3,isite)
         dsqr=x**2+y**2+z**2
c        --- reciprocals of d
         d1=one/sqrt(dsqr)
         d3=d1/dsqr
c        --- get contribution of site
         vcalcr=vcalcr+
     $         zm(isite)*(zl(1,isite)*x+zl(2,isite)*y+zl(3,isite)*z)*d3
   10 continue
c
c
      return
      end
