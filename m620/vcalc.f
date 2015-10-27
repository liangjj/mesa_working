*deck @(#)vcalc.f	5.1  11/6/94
      function vcalc(nsite,zm,site,vo)
c
      implicit none
      real*8 vcalc
      integer nsite
      real*8 zm(nsite,10),site(3,nsite),vo(4)
c
      integer isite
      real*8 x,y,z,dsqr,d1,d3,d5
      real*8 zero,one,half,three
c
      parameter (zero=0.0d+00,one=1.0d+00,half=0.5d+00,three=3.0d+00)
c
c
      vcalc=zero
      do 10 isite=1,nsite
         x=vo(1)-site(1,isite)
         y=vo(2)-site(2,isite)
         z=vo(3)-site(3,isite)
         dsqr=x**2+y**2+z**2
c        --- reciprocals of d
         d1=one/sqrt(dsqr)
         d3=d1/dsqr
         d5=d3/dsqr
c        --- get contribution of site
         vcalc=vcalc+zm(isite,1)*d1
     $           +zm(isite,2)*x*d3
     $           +zm(isite,3)*y*d3
     $           +zm(isite,4)*z*d3
     $           +zm(isite,5)*half*x*x*d5
     $           +zm(isite,6)*half*y*y*d5
     $           +zm(isite,7)*half*z*z*d5
     $           +zm(isite,8)*three*x*y*d5
     $           +zm(isite,9)*three*x*z*d5
     $           +zm(isite,10)*three*y*z*d5
   10 continue
c
c
      return
      end
