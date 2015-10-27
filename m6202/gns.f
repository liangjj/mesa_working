*deck gns
c***begin prologue     gns
c***date written       940811   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           gaussian, link m6202
c***author             schneider, barry (nsf)
c***source             m6202
c***purpose            gaussian on a grid 
c***references         none
c***routines called
c***end prologue        gns
      subroutine gns(f,nx,ny,nz,alpha,cen,r,ctheta,cphi,sphi,
     1               stheta,nshell,nrtot,nr,nthet,nphi)
      implicit integer (a-z)
      real*8 f, alpha, cen, r, ctheta, stheta, cphi, sphi
      real*8 x, xc, y, yc, z, zc
      real*8 zz, dsq, zc2
      dimension f(nrtot,nthet,nphi), r(nrtot), ctheta(nthet)
      dimension stheta(nthet), sphi(nphi), cphi(nphi)
      dimension cen(3), nr(nshell)
      common/io/ inp,iout
      do 10 thept=1,nthet
         stheta(thept)=sqrt(1.d0-ctheta(thept)*ctheta(thept))
   10 continue
      rcount=0
      do 20 ns=1,nshell
         do 30 rpt=1,nr(ns)
            rcount=rcount+1
            do 40 thept=1,nthet
               z=r(rcount)*ctheta(thept)
               zc=z-cen(3)
               zc2=zc*zc
               zc=zc**nz
               zz=r(rcount)*stheta(thept)
               do 50 phipt=1,nphi
                  x=zz*cphi(phipt)
                  xc=(x-cen(1))**nx
                  y=zz*sphi(phipt)
                  yc=(y-cen(2))**ny
                  dsq=(x-cen(1))*(x-cen(1)) +
     1                (y-cen(2))*(y-cen(2)) + 
     2                         zc2
                  f(rcount,thept,phipt) = xc*yc*zc*exp(-alpha*dsq)
   50          continue                                           
   40       continue
   30    continue
   20 continue
      return
      end






