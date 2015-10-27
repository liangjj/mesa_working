*deck dgsn
c***begin prologue     dgsn
c***date written       940811   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           yukawa, link m6201
c***author             schneider, barry (nsf)
c***source             m6202
c***purpose            derivative of gaussian on a grid with respect
c***                   to radial co-ordinate.
c***references         none
c                      storage of fun is packed. when a lebedev quadrature
c                      is being used nthet=nphi=nang and the values of all
c                      angular quantities are always used together.
c***routines called
c***end prologue        dgsn
      subroutine dgsn(dfun,nx,ny,nz,alpha,cen,r,ctheta,cphi,sphi,
     1                stheta,nshell,nrtot,nr,nthet,nphi,
     2                nang,nonsep)
      implicit integer (a-z)
      real*8 dfun, alpha, cen, r, rad, ctheta, stheta, cphi, sphi
      real*8 x, xc, y, yc, z, zc
      real*8 xcp, xcm, ycp,ycm, zcp, zcm
      real*8 xovr, yovr, zovr
      real*8 zz, dsq, zc2
      logical nonsep
      dimension dfun(nrtot,*), r(nrtot), ctheta(nthet)
      dimension stheta(nthet), sphi(nphi), cphi(nphi)
      dimension cen(3), nr(nshell)
      common/io/ inp,iout
      nprd=nthet*nphi
      if (nonsep) then
          nprd=nang
      endif     
      do 10 thept=1,nthet
         stheta(thept)=sqrt(1.d0-ctheta(thept)*ctheta(thept))
   10 continue
      rcount=0
      if (nonsep) then      
          do 20 ns=1,nshell
             do 30 rpt=1,nr(ns)
                rcount=rcount+1       
                do 40 ang=1,nang
                   z=r(rcount)*ctheta(ang)
                   zc=(z-cen(3))**nz
                   zcp=zc*(z-cen(3))
                   zcm=zc/(z-cen(3))
                   zz=r(rcount)*stheta(ang)
                   x=zz*cphi(ang)
                   xc=(x-cen(1))**nx
                   xcp=(x-cen(1))*xc
                   xcm=xc/(x-cen(1))
                   y=zz*sphi(ang)      
                   yc=(y-cen(2))**ny
                   ycp=(y-cen(2))*yc
                   ycm=yc/(y-cen(2))
                   rad=sqrt(x*x+y*y+z*z)
                   zovr=z/rad
                   xovr=x/rad
                   yovr=y/rad
                   dsq=(x-cen(1))*(x-cen(1)) +
     1                 (y-cen(2))*(y-cen(2)) + 
     2                 (z-cen(3))*(z-cen(3))
                   dfun(rcount,ang) = ( -2.d0*alpha*( xovr*xcp*yc*zc +
     1                                                yovr*xc*ycp*zc +
     2                                                zovr*xc*yc*zcp ) +
     3                                                nx*xcm*yc*zc     +
     4                                                ny*xc*ycm*zc     +
     5                                                nz*xc*yc*zcm )*
     6                                                exp(-alpha*dsq)
   40           continue                                           
   30        continue
   20     continue
      else
          do 60 ns=1,nshell
             do 70 rpt=1,nr(ns)
                rcount=rcount+1
                ang=0       
                do 80 thept=1,nthet
                   z=r(rcount)*ctheta(thept)
                   zc=z-cen(3)
                   zc2=zc*zc
                   zc=zc**nz
                   zcp=zc*(z-cen(3))
                   zcm=zc/(z-cen(3))
                   zovr=z/rad
                   zz=r(rcount)*stheta(thept)
                   do 90 phipt=1,nphi
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
     1                    (y-cen(2))*(y-cen(2)) + 
     2                           zc2
                      ang=ang+1      
                      dfun(rcount,ang) = ( -2.d0*alpha*
     1                                      ( xovr*xcp*yc*zc +
     2                                        yovr*xc*ycp*zc +
     3                                        zovr*xc*yc*zcp ) +
     4                                        nx*xcm*yc*zc + 
     5                                        ny*xc*ycm*zc +
     6                                        nz*xc*yc*zcm )*
     5                                        exp(-alpha*dsq)
   90              continue                                           
   80           continue
   70        continue
   60     continue
      endif
      return
      end






