*deck @(#)newf2ylm.f	5.1  11/28/95
      subroutine newf2ylm(nrblk,angmx,nlmx,rblkord,rblkbgn,nrblk,
     $     rblksiz,nlm  )
c***begin prologue     newf2ylm.f
c***date written       951121      (yymmdd)  
c***revision date      11/28/95
c
c***keywords           
c***author             martin, richard(lanl); russo, thomas(lanl)
c***source             @(#)newf2ylm.f	5.1 11/28/95
c***purpose            performs decomposition of a function into 
c                      spherical harmonic components on general or pruned
c                      grids.
c***description        
c
C    for each block of radial shells on this atom with like angular grids,
c    read in the appropriate set of ylms, fit each radial shell in this
c    block to ylms up to the appropriate order.
c    Some rholms will wind up being zero at spots because of the differing
c    orders.
c    This is similar to the "ftoylm" subroutine, but that routine assumed 
c    that all radii had same angular quadrature.  This one uses some
c    additional data structures to represent the grid so that each radii
c    can have a different angular grid.
c
c***references
c
c***routines called
c
c***end prologue       newf2ylm.f
      implicit none
c --- input variables
c     nrblk: number of radial shell blocks of same lebedev type
c     angmx: number of angular points for largest grid angular in this atom
c     nlmx: largest number of ylms on any grid in this atom
      integer nrblk,angmx,nlmx
c --- input arrays (unmodified)
c     rblkord: lebedev orders of radial blocks
c     rblkbgn: starting radial point for given block
c     rblksiz: number of shells in given block
c     nlm: number of ylms to do for given block (depends on lebord)
      integer rblkord(nrblk),rblkbgn(nrblk),rblksiz(nrblk)
      integer nlm(nrblk)
c --- input arrays (tweaked)
c --- output variables
c --- output arrays
c --- scratch arrays
c --- typed functions
      real*8 sdot
c --- local variables
      integer irshlb,lebord,rmin,rmax,ir,ilm

      character*32 grdnm
 1100 format('atom ',i3.3,' angwt ',i3.3)
 1200 format('atom ',i3.3,' ylm ',i3.3)

      call rzero(flm,nlmx*nr)
      do 100 irshlb=1,nrblk
         lebord=rblkord(irshlb)
         rmin=rblkbgn(irshlb)
         rmax=rblkbgn(irshlb)+rblksiz(irshlb)-1
         write(unit=grdnm,fmt=1100),iatom,irshlb
c load in wts for this lebord into angwts, size into angsiz
         call iosys('read real "'//grdnm//'" from '//grdfil,
     $        -1,angwt,0,' ')
c load in the ylms for this lebord (generate?  read in?) into ylm(angsiz,nlm)
         write(unit=grdnm,fmt=1200),iatom,irshlb
         call iosys('read real "'//grdnm//'" from '//grdfil,
     $        -1,ylm,0,' ')
         do 200 ir=rmin,rmax
c           multiply function on radial shell by angular weights, then
c           multiply by ylms and sum over grid, thereby forming
c                  flm(r,lm)=integral(Y(omega,lm)*f(r,omega))d(omega)
            call vmul(scr,angwt,f(1,ir),angsiz)
            do 300 ilm=1,nlm(irshlb)
               flm(ir,ilm)=sdot(angsiz,scr,1,ylm(1,ilm),1)
 300        continue 
 200     continue 
 100  continue 

      return
      end
