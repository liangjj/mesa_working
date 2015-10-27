*deck @(#)ftoylm.f	1.1  4/25/95
      subroutine ftoylm(f,flm,ylm,nr,nomega,nlm,wt,scr,ylmerr)
c***begin prologue     ftoylm.f
c***date written       940304  
c***revision date      11/6/94      
c
c***keywords           
c***author             martin, richard(lanl) 
c***source             @(#)ftoylm.f	1.1   4/25/95
c***purpose            performs decomposition of a function into 
c                      spherical harmonic components.
c***description
c   f       ...  (nomega,nr)  the function to be decomposed
c   flm     ...  (nr,nlm)   decomposition into lm components
c   ylm     ...  (nomega,nlm) real spherical harmonics on the unit sphere.
c   nr      ...  the number of radial points.
c   nomega    ...  the number of angular points.
c   nlm     ...  the total number of ylm functions.
c   wt      ...  the weights(nomega).
c   scr     ...  scratch arrray(nomega,nr)
c   ylmerr  ...  this is used as a measure of "zero". it is the largest
c                deviation from orthonormality for the ylm.
c    
c***references
c
c***routines called
c
c***end prologue       ftoylm.f
      implicit none
c     --- input variables -----
      integer nlm,nr,nomega
      real*8 ylmerr
c     --- input arrays (unmodified) ---
      real*8 f(nomega,nr),wt(nomega),ylm(nomega,nlm)
c     --- input arrays (scratch) ---
      real*8 scr(nomega,nr)
c     --- output arrays ---
      real*8 flm(nr,nlm)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer shell
      real*8 zero,one
      parameter (zero=0.0d+00,one=1.0d+00)
c
c     --- multiply the function by the weights and integrate
      do 10 shell=1,nr
         call vmul(scr(1,shell),wt,f(1,shell),nomega)
   10 continue
      call ebtc(flm,scr,ylm,nr,nomega,nlm)
c
c     --- set all components of flm which are smaller than ylmerr to zero.
      call vclean(flm,ylmerr,nr*nlm)
c
c
c     --- note on use:
c         to use variable angular meshes for the radial coordinate,
c         one will have to have a loop in the calling routine
c         pass it only the points which are associated with the same level
c         of angular quadrature. 
c         do shell=1,radshls
c            nr = number of shells with same lebedev order.
c            call lebedev to get wts. (if the original weights are used,
c                then function will have radial weights included.
c            call ylm to get functions.
c            call decomposition scheme.
c         enddo
c
c
      return
      end
