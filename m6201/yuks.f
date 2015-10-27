*deck yuks
c***begin prologue     yuks
c***date written       930524   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           yukawa, link m6201
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            calculate multicenter yukawa potential
c***references         none
c                      storage of fun is packed. when a lebedev quadrature
c                      is being used nthet=nphi=nang and the values of all
c                      angular quantities are always used together.
c***routines called
c***end prologue        yuks
      subroutine yuks(f,alpha,cen,r,ctheta,cphi,sphi,stheta,
     1                ncen,nshell,nrtot,nr,nthet,nphi)
      implicit integer (a-z)
      real*8 f, funs, alpha, cen, r, ctheta, stheta, cphi, sphi
      real*8 x, y, z, zz, d, pi, fourpi
      dimension f(nrtot,nthet,nphi)
      dimension ctheta(nthet), stheta(nthet), sphi(nphi), cphi(nphi)
      dimension cen(3,ncen), alpha(ncen), nr(nshell), r(nrtot)
      common/io/ inp,iout
      data pi / 3.14159265358979323846d+00 /
      fourpi=4.d0*pi
      do 10 thept=1,nthet
         stheta(thept)=sqrt(1.d0-ctheta(thept)*ctheta(thept))
   10 continue
      rcount=0
      call rzero(funs,nrtot*nthet*nphi)
      do 20 ns=1,nshell
         do 30 rpt=1,nr(ns)
            rcount=rcount+1
            do 40 thept=1,nthet
               z=r(rcount)*ctheta(thept)
               zz=r(rcount)*stheta(thept)
               do 50 phipt=1,nphi
                  x=zz*cphi(phipt)
                  y=zz*sphi(phipt)
                  do 60 icen=1,ncen
                     d =sqrt( (x-cen(1,icen))*(x-cen(1,icen)) +
     1                        (y-cen(2,icen))*(y-cen(2,icen)) + 
     2                        (z-cen(3,icen))*(z-cen(3,icen)) )
c                     f(rcount,thept,phipt) = 
c     1                         f(rcount,thept,phipt) + 
c     2                           exp(-alpha(icen)*d)/d
                     f(rcount,thept,phipt) = 
     1                         f(rcount,thept,phipt) +
     2                           d*exp(-2.d0*alpha(icen)*d)
   60             continue
c                      f(rcount,thept,phipt)=z*exp(-(x*x+y*y+z*z))
   50          continue                                           
   40       continue
   30    continue
   20 continue
      call sscal(rcount*nthet*nphi,1.d0/pi,f,1)
      call sscal(rcount*nthet*nphi,-fourpi,f,1)
      return
      end
