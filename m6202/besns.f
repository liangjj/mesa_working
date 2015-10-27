*deck besns.f
c***begin prologue     besns
c***date written       940811   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           yukawa, link m6201
c***author             schneider, barry (nsf)
c***source             m6202
c***purpose            ricatti-bessel on a grid 
c***references         none
c***routines called
c***end prologue        besns
      subroutine besns(f,k,l,j,xval,cen,r,ctheta,cphi,sphi,
     1                 stheta,nshell,nrtot,nr,nthet,nphi,
     2                 nang,ltop)
      implicit integer (a-z)
      parameter (acc=30)
      real*8 f, k, cen, r, ctheta, stheta, cphi, sphi
      real*8 x, xc, y, yc, z, zc
      real*8 zz, dsq, zc2
      real*8 j, jp, yy, yp, xval 
      dimension f(nrtot,nang), r(nrtot), ctheta(nthet)
      dimension stheta(nthet), sphi(nphi), cphi(nphi)
      dimension cen(3), nr(nshell), xval(nang), j(nang,0:ltop)
      common/io/ inp,iout
      do 10 thept=1,nthet
         stheta(thept)=sqrt(1.d0-ctheta(thept)*ctheta(thept))
   10 continue
      rcount=0
      do 20 ns=1,nshell
         do 30 rpt=1,nr(ns)
            rcount=rcount+1       
            do 40 ang=1,nang
               z=r(rcount)*ctheta(ang)
               zc=(z-cen(3))**nz
               zz=r(rcount)*stheta(ang)
               x=zz*cphi(ang)
               xc=(x-cen(1))**nx
               y=zz*sphi(ang)      
               yc=(y-cen(2))**ny
               dsq=(x-cen(1))*(x-cen(1)) +
     1             (y-cen(2))*(y-cen(2)) + 
     2             (z-cen(3))*(z-cen(3))
               xval(ang) = k*sqrt(dsq)
   40       continue                                           
            call rcbssl(xval,j,jp,yy,yp,nang,l,ltop,
     1                  'no derivatives','regular',.false.)
            do 50 ang=1,nang
               f(rcount,ang)=j(ang,l)
   50       continue
   30    continue
   20 continue
      return
      end








