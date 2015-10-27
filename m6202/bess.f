*deck bess.f
c***begin prologue     bess
c***date written       940811   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           yukawa, link m6201
c***author             schneider, barry (nsf)
c***source             m6202
c***purpose            ricatti-bessel on a grid 
c***references         none
c***routines called
c***end prologue        bess
      subroutine bes(f,k,l,j,xval,cen,r,ctheta,cphi,sphi,
     1               stheta,nshell,nrtot,nr,nthet,nphi,
     2               ltop)
      implicit integer (a-z)
      parameter (acc=30)
      real*8 f, k, cen, r, ctheta, stheta, cphi, sphi
      real*8 x, xc, y, yc, z, zc
      real*8 zz, dsq, zc2
      real*8 j, jp, yy, yp, xval 
      dimension f(nrtot,nthet,nphi), r(nrtot), ctheta(nthet)
      dimension stheta(nthet), sphi(nphi), cphi(nphi)
      dimension cen(3), nr(nshell), xval(nthet*nphi)
      dimension j(nthet*nphi,0:ltop)
      common/io/ inp,iout
      do 10 thept=1,nthet
         stheta(thept)=sqrt(1.d0-ctheta(thept)*ctheta(thept))
   10 continue
      rcount=0
      do 20 ns=1,nshell
         do 30 rpt=1,nr(ns)
            rcount=rcount+1
            angl=0
            do 40 thept=1,nthet
               z=r(rcount)*ctheta(thept)
               zc=z-cen(3)
               zc2=zc*zc
               zc=zc**nz
               zz=r(rcount)*stheta(thept)
               do 50 phipt=1,nphi
                  angl=angl+1
                  x=zz*cphi(phipt)
                  xc=(x-cen(1))**nx
                  y=zz*sphi(phipt)
                  yc=(y-cen(2))**ny
                  dsq=(x-cen(1))*(x-cen(1)) +
     1                (y-cen(2))*(y-cen(2)) + 
     2                       zc2
                  xval(angl) = k*sqrt(dsq)
   50          continue                                           
   40       continue
            call rcbssl(xval,j,jp,yy,yp,angl,l,ltop,
     1                  'no derivatives','regular',.false.)
            angl=0 
            do 60 thept=1,nthet
               do 70 phipt=1,nphi
                  angl=angl+1
                  f(rcount,thept,phipt)=j(angl,l)
   70          continue
   60       continue
   30    continue         
   20 continue
      return
      end








