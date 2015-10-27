*deck wdcalc.f
c***begin prologue     wdcalc
c***date written       941019   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6236, link 6236, propagators
c***author             schneider, b. i.(nsf)
c***source             m6236
c***purpose            calculation of propagators
c***description        
c***                   
c***                   
c***references         papers by light-walker j. chem. phys.
c
c***routines called    iosys, util and mdutil
c***end prologue       m6236
c                    
      subroutine wdcalc (en,du,d,e,nc)
c     calculate the sector propagator 
c
c                     du = left - right r value
c                     ( negative of step size moving to larger r values )
c
      implicit integer(a-z)
      real*8 en, du, d, e, ek, hek, s, sh, ch
      common/io/inp,iout
      dimension d(nc,4), en(nc)
      do 10 i=1,nc
         ek=e-en(i)
         if (ek.ge.0.d0) then
             ek=sqrt(ek)
             hek=du*ek
             s=-1.d0
             sh=sin(hek)
             ch=cos(hek)
         else   
             s=1.d0
             ek=sqrt(-ek)
             hek=du*ek
             sh=sinh(hek)
             ch=cosh(hek)
         endif
         d(i,1)=ch
         d(i,2)=sh/ek
         d(i,3)=s*sh*ek
         d(i,4)=ch
 10   continue    
      return
      end
