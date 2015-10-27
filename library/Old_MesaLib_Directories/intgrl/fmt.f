*deck @(#)fmt.f	5.1  11/6/94
      subroutine fmt(prmint,xyz,npint,lenblk,imax,jmax,t1,mini,maxi,
     $               minj,maxj,nx,ny,nz)
c***begin prologue     fmt.f
c***date written       840723  
c***revision date      11/6/94      
c
c***keywords           
c***author             saxe, paul(lanl)
c***source             @(#)fmt.f	5.1   11/6/94
c***purpose            
c***description
c      module to assemble the two-dimensional integrals in xyz into
c      primitive kinetic-energy integrals in prmint.
c     
c***references
c
c***routines called
c
c***end prologue       fmt.f
      implicit none
c     --- input variables -----
      integer npint,lenblk,imax,jmax,mini,maxi,minj,maxj
c     --- input arrays (unmodified) ---
      integer nx(*),ny(*),nz(*)
      real*8 xyz(npint,0:imax,0:jmax,3,2)
c     --- input arrays (scratch) ---
      real*8 t1(npint)
c     --- output arrays ---
      real*8 prmint(npint,lenblk)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer intgrl,i,j,ix,iy,iz,jx,jy,jz
c
c
      intgrl=0
c
      do 2 i=mini,maxi
         ix=nx(i)
         iy=ny(i)
         iz=nz(i)
         do 1 j=minj,maxj
            jx=nx(j)
            jy=ny(j)
            jz=nz(j)
            intgrl=intgrl+1
c
            call vmul(prmint(1,intgrl),xyz(1,ix,jx,1,2),
     $                                 xyz(1,iy,jy,2,1),npint)
            call vmul(t1,xyz(1,ix,jx,1,1),xyz(1,iy,jy,2,2),npint)
            call vadd(prmint(1,intgrl),prmint(1,intgrl),t1,npint)
            call vmul(prmint(1,intgrl),prmint(1,intgrl),
     $                                 xyz(1,iz,jz,3,1),npint)
            call vmul(t1,xyz(1,ix,jx,1,1),xyz(1,iy,jy,2,1),npint)
            call vmul(t1,t1,xyz(1,iz,jz,3,2),npint)
            call vadd(prmint(1,intgrl),prmint(1,intgrl),t1,npint)
    1    continue
    2 continue
c
c
      return
      end
