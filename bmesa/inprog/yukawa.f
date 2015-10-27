*deck @(#)yukawa.f	1.1 9/9/91
c***begin prologue     yukawa
c***date written       930524   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           yukawa, link m6201
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            calculate multicenter yukawa potential
c***references         none
c                      storage of fun is packed.
c***routines called
c***end prologue        yukawa
      subroutine yukawa(fun,alpha,cen,r,theta,cphi,sphi,ncen,nr,
     1                  nthet,nphi)
      implicit integer (a-z)
      real*8 fun, alpha, cen, r, theta, cphi, sphi, snthe, temp
      real*8 x, y, z, d
      dimension fun(nr,nthet,nphi), r(nr), theta(nthet), cphi(nphi)
      dimension sphi(nphi), cen(3,ncen), alpha(ncen)
      common/io/ inp,iout
      call rzero(fun,nr*nthet*nphi)
      do 10 thept=1,nthet
         snthe=sqrt(1.d0-theta(thept)*theta(thept))
         do 20 rpt=1,nr
            temp=r(rpt)*snthe    
            z=r(rpt)*theta(thept)
            do 30 phipt=1,nphi
               x=temp*cphi(phipt)
               y=temp*sphi(phipt)      
               do 40 icen=1,ncen
                  d =sqrt( (x-cen(1,icen))*(x-cen(1,icen))+
     1                     (y-cen(2,icen))*(y-cen(2,icen)) + 
     2                     (z-cen(3,icen))*(z-cen(3,icen)) )
                  fun(rpt,thept,phipt) = fun(rpt,thept,phipt) + 
     1                                   exp(-alpha(icen)*d)/d
  40           continue
  30        continue
  20     continue                                           
  10  continue
      return
      end
