*deck dgns
c***begin prologue     dgns
c***date written       940811   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           gaussian, link m6201
c***author             schneider, barry (nsf)
c***source             m6202
c***purpose            derivative of gaussian on a grid with respect
c***                   to radial co-ordinate.
c***references         none
c***routines called
c***end prologue        dgns
      subroutine dgns(df,nx,ny,nz,alpha,cen,r,ctheta,cphi,sphi,
     1                stheta,nshell,nrtot,nr,nthet,nphi)
      implicit integer (a-z)
      real*8 df, alpha, cen, r, rad, ctheta, stheta, cphi, sphi
      real*8 x, xc, y, yc, z, zc
      real*8 xcp, xcm, ycp,ycm, zcp, zcm
      real*8 xovr, yovr, zovr
      real*8 zz, dsq, zc2
      dimension df(nrtot,nthet,nphi), r(nrtot), ctheta(nthet)
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
               zcp=zc*(z-cen(3))
               zcm=zc/(z-cen(3))
               zovr=z/rad
               zz=r(rcount)*stheta(thept)
               do 50 phipt=1,nphi
                  x=zz*cphi(phipt)
                  xc=(x-cen(1))**nx
                  xcp=(x-cen(1))*xc
                  xcm=xc/(x-cen(1))
                  y=zz*sphi(phipt)
                  yc=(y-cen(2))**ny
                  ycp=(y-cen(2))*yc
                  ycm=yc/(y-cen(2))
                  rad=sqrt(x*x+y*y+z*z)
                  xovr=x/rad
                  yovr=y/rad
                  dsq=(x-cen(1))*(x-cen(1)) +
     1                (y-cen(2))*(y-cen(2)) + 
     2                       zc2
                  df(rcount,thept,phipt) = ( -2.d0*alpha*
     1                                     ( xovr*xcp*yc*zc +
     2                                       yovr*xc*ycp*zc +
     3                                       zovr*xc*yc*zcp ) +
     4                                       nx*xcm*yc*zc + 
     5                                       ny*xc*ycm*zc +
     6                                       nz*xc*yc*zcm )*
     5                                       exp(-alpha*dsq)
   50          continue                                           
   40       continue
   30    continue
   20 continue
      return
      end






