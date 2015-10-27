c $Header$
*deck pertsl.f
c***begin prologue     pertsl
c***date written       930201   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           pertsl, link 6201, wavefunction
c***author             schneider, barry (nsf)
c***source             m6203
c***purpose            solve inhomogeneous radial schroedinger equation
c*** 
c***description        the first order inhomogeneous radial schroedinger
c***                   equation for an exponentially damped bessel function
c***                   is solved as an integral equation. the solution is
c***                   used as a basis function in the variational
c***                   calculation.
c***                    
c***references         
c***routines called
c***end prologue       pertsl
      subroutine pertsl (psln,greg,gireg,pt,wt,rhsfn,nr,nl,prnt)
      implicit integer(a-z)
      real*8 psln, greg, gireg, pt, wt, rhsfn, sumf, sumb
      character*4 fcall
      logical prnt
      common /io/ inp, iout
      dimension psln(nr,nl), greg(nr,nl), gireg(nr,nl)
      dimension pt(nr), wt(nr), rhsfn(nr,nl)
*
*
*          the radial functions are defined as solutions
*          of the equation;
*            g'' + ( 2*abs(e) - l*(l+1)/(r+r) ) g  = -2.d0*exp(-r)*greg(kr)
*                -
*
*
      do 10 ipt=1,nr
         rhsfn(ipt,1)=-2.d0*exp(-pt(ipt))
   10 continue
      do 20 l=2,nl
         call copy(rhsfn(1,1),rhsfn(1,l),nr)
   20 continue 
      do 20 l=1,nl
         call vmul(rhsfn(1,l),rhsfn(1,l),greg(1,l),nr)
         sumf=0.d0
         do 30 ipt=1,nr
            sumf=sumf+greg(ipt,l)*wt(ipt)*rhsfn(ipt,l)
            psln(ipt,l)=gireg(ipt,l)*sumf
   30    continue
         sumb=0.d0
         do 40 ipt=nr-1,1,-1
            ip=ipt+1
            sumb=sumb+gireg(ip,l)*wt(ip)*rhsfn(ip,l)
            psln(ipt,l)=psln(ipt,l)+greg(ipt,l)*sumb
   40    continue
   20 continue
      return
      end


