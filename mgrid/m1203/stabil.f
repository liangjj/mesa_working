*deck stabil
c***begin prologue     stabil
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            stability analysis for three dimensional partial
c***                   differential equation.
c***description 
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       stabil
      subroutine stabil(stp,k)
c
      implicit integer (a-z)
      real*8 stp, k, fac  
      fac=6.d0-k*k*stp*stp
      if (fac.ge.5.d0) then
          if (fac.le.6.d0) then
              write(iout,1) fac
          else
              write(iout,2) fac
          endif                  
      endif
      return
    1 format(/,1x,'outside stable region.  fac = ',e15.8)
    2 format(/,1x,'inside stable region.  fac = ',e15.8)     
      end





