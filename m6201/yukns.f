*deck yukns
c***begin prologue     yukns
c***date written       930524   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           yukns, link m6201
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            calculate multicenter yukawa potential
c***references         none
c***routines called
c***end prologue        yukns
      subroutine yukns(f,alpha,cen,r,ctheta,cphi,sphi,stheta,
     1                 ncen,nshell,nrtot,nr,nang)
      implicit integer (a-z)
      real*8 f, alpha, cen, r, ctheta, stheta, cphi, sphi
      real*8 x, y, z, zz, d, pi, fourpi
      dimension f(nrtot,nang)
      dimension ctheta(nang), stheta(nang), sphi(nang), cphi(nang)
      dimension cen(3,ncen), alpha(ncen), r(nrtot), nr(nshell)
      common/io/ inp,iout
      data pi / 3.14159265358979323846d+00 /
      fourpi=4.d0*pi
      do 10 thept=1,nang
         stheta(thept)=sqrt(1.d0-ctheta(thept)*ctheta(thept))
   10 continue
      rcount=0
      call rzero(f,nrtot*nang)      
      do 20 ns=1,nshell
         do 30 rpt=1,nr(ns)
            rcount=rcount+1       
            do 40 ang=1,nang
               z=r(rcount)*ctheta(ang)
               zz=r(rcount)*stheta(ang)
               x=zz*cphi(ang)
               y=zz*sphi(ang)      
               do 50 icen=1,ncen
                  d =sqrt( (x-cen(1,icen))*(x-cen(1,icen)) +
     1                     (y-cen(2,icen))*(y-cen(2,icen)) + 
     2                     (z-cen(3,icen))*(z-cen(3,icen)) )
                  f(rcount,ang)=f(rcount,ang) +
     1                            d*exp(-2.d0*alpha(icen)*d)
c                  f(rcount,ang) = f(rcount,ang) + 
c     1                              exp(-alpha(icen)*d)/d
c                  f(rcount,ang)=z*exp(-(x*x+y*y+z*z))
   50          continue
   40       continue                                           
   30    continue
   20 continue
      call sscal(rcount*nang,1.d0/pi,f,1)
      call sscal(rcount*nang,-fourpi,f,1)
      return
      end
