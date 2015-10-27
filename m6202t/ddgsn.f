*deck ddgsn
c***begin prologue     ddgsn
c***date written       940811   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           yukawa, link m6201
c***author             schneider, barry (nsf)
c***source             m6202
c***purpose            laplacian of gaussian on a grid
c***references         none
c                      storage of fun is packed. when a lebedev quadrature
c                      is being used nthet=nphi=nang and the values of all
c                      angular quantities are always used together.
c***routines called
c***end prologue        ddgsn
      subroutine ddgsn(ddfun,nx,ny,nz,alpha,cen,r,ctheta,cphi,sphi,
     1                stheta,nshell,nrtot,nr,nthet,nphi,
     2                nang,nonsep)
      implicit integer (a-z)
      real*8 ddfun, alpha, cen, r, ctheta, stheta, cphi, sphi
      real*8 x, xc, y, yc, z, zc
      real*8 xcpp, xcmm, ycpp, ycmm, zcpp, zcmm
      real*8 zz, dsq, zc2
      logical nonsep
      dimension ddfun(nrtot,*), r(nrtot), ctheta(nthet)
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
                   zcpp=zc*(z-cen(3))*(z-cen(3))
                   zcmm=zc/( (z-cen(3))*(z-cen(3)) )
                   zz=r(rcount)*stheta(ang)
                   x=zz*cphi(ang)
                   xc=(x-cen(1))**nx
                   xcpp=(x-cen(1))*(x-cen(1))*xc
                   xcmm=xc/( (x-cen(1))*(x-cen(1)) )
                   y=zz*sphi(ang)      
                   yc=(y-cen(2))**ny
                   ycpp=(y-cen(2))*(y-cen(2))*yc
                   ycmm=yc/( (y-cen(2))*(y-cen(2)) )
                   dsq=(x-cen(1))*(x-cen(1)) +
     1                 (y-cen(2))*(y-cen(2)) + 
     2                 (z-cen(3))*(z-cen(3))
                   ddfun(rcount,ang) = ( nx*(nx-1)*xcmm*yc*zc     +
     1                                   ny*(ny-1)*xc*ycmm*zc     +
     2                                   nz*(nz-1)*xc*yc*zcmm     -
     3                                   2.d0*alpha*(nx+nx+ny+ny+
     4                                   nz+nz+3)*xc*yc*zc        +
     5                                   4.d0*alpha*alpha*
     6                                   (xcpp*yc*zc+xc*ycpp*zc +
     7                                    xc*yc*zcpp ) )*exp(-alpha*dsq)
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
                   zcpp=zc*(z-cen(3))*(z-cen(3))
                   zcmm=zc/( (z-cen(3))*(z-cen(3)) )
                   zz=r(rcount)*stheta(thept)
                   do 90 phipt=1,nphi
                      x=zz*cphi(phipt)
                      xc=(x-cen(1))**nx
                      xcpp=(x-cen(1))*(x-cen(1))*xc
                      xcmm=xc/( (x-cen(1))*(x-cen(1)) )
                      y=zz*sphi(phipt)
                      yc=(y-cen(2))**ny
                      ycpp=(y-cen(2))*(y-cen(2))*yc
                      ycmm=yc/( (y-cen(2))*(y-cen(2)) )
                      dsq=(x-cen(1))*(x-cen(1)) +
     1                    (y-cen(2))*(y-cen(2)) + 
     2                           zc2
                      ang=ang+1
                      ddfun(rcount,ang) = ( nx*(nx-1)*xcmm*yc*zc     +
     1                                      ny*(ny-1)*xc*ycmm*zc     +
     2                                      nz*(nz-1)*xc*yc*zcmm     -
     3                                      2.d0*alpha*(nx+nx+ny+ny+
     4                                      nz+nz+3)*xc*yc*zc        +
     5                                      4.d0*alpha*alpha*
     6                                      (xcpp*yc*zc+xc*ycpp*zc +
     7                                       xc*yc*zcpp) )*
     8                                                 exp(-alpha*dsq)      
   90              continue                                           
   80           continue
   70        continue
   60     continue
      endif
      return
      end






