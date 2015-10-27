*deck yukawa
c***begin prologue     yukawa
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
c***end prologue        yukawa
      subroutine yukawa(fun,alpha,cen,r,ctheta,cphi,sphi,stheta,ncen,
     1                  nshell,nrtot,nr,nthet,nphi,nang,nonsep)
      implicit integer (a-z)
      real*8 fun, alpha, cen, r, ctheta, stheta, cphi, sphi
      real*8 x, y, z, zz, d, pi, fourpi
      logical nonsep
      dimension fun(nrtot,*), r(nrtot), ctheta(nthet)
      dimension stheta(nthet), sphi(nphi), cphi(nphi), cen(3,ncen)
      dimension alpha(ncen), nr(nshell)
      common/io/ inp,iout
      data pi / 3.14159265358979323846d+00 /
      fourpi=4.d0*pi
      nprd=nthet*nphi
      if (nonsep) then
          nprd=nang
      endif     
      call rzero(fun,nrtot*nprd)
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
                   zz=r(rcount)*stheta(ang)
                   x=zz*cphi(ang)
                   y=zz*sphi(ang)      
                   do 50 icen=1,ncen
                      d =sqrt( (x-cen(1,icen))*(x-cen(1,icen)) +
     1                         (y-cen(2,icen))*(y-cen(2,icen)) + 
     2                         (z-cen(3,icen))*(z-cen(3,icen)) )
                      fun(rcount,ang)=fun(rcount,ang) +
     1                                 d*exp(-2.d0*alpha(icen)*d)
c                      fun(rcount,ang) = fun(rcount,ang) + 
c     1                                       exp(-alpha(icen)*d)/d
   50              continue

   40           continue                                           
   30        continue
   20     continue
      else
          do 60 ns=1,nshell
             do 70 rpt=1,nr(ns)
                rcount=rcount+1
                angcnt=0       
                do 80 thept=1,nthet
                   z=r(rcount)*ctheta(thept)
                   zz=r(rcount)*stheta(thept)
                   do 90 phipt=1,nphi
                      x=zz*cphi(phipt)
                      y=zz*sphi(phipt)
                      angcnt=angcnt+1      
                      do 100 icen=1,ncen
                         d =sqrt( (x-cen(1,icen))*(x-cen(1,icen)) +
     1                            (y-cen(2,icen))*(y-cen(2,icen)) + 
     2                            (z-cen(3,icen))*(z-cen(3,icen)) )
c                         fun(rcount,angcnt) = fun(rcount,angcnt) + 
c     1                                        exp(-alpha(icen)*d)/d
                         fun(rcount,angcnt)=fun(rcount,angcnt) +
     1                                      d*exp(-2.d0*alpha(icen)*d)
  100                 continue
   90              continue                                           
   80           continue
   70        continue
   60     continue
      endif
      call sscal(nprd*rcount,1.d0/pi,fun,1)
      call sscal(nprd*rcount,-fourpi,fun,1)          
      return
      end
