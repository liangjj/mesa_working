*deck @(#)fmdel.f	5.1  11/6/94
      subroutine fmdel(prmint,xyz,npint,lenblk,imax,jmax,mini,maxi,
     $                  minj,maxj,nx,ny,nz,coord)
c***begin prologue     fmdel.f
c***date written       840723  
c***revision date      11/6/94      
c   september 8,1986   rlm at lanl
c      modifying fmonel to handle momentum operator.
c
c***keywords           
c***author             saxe, paul and martin, richard(lanl) 
c***source             @(#)fmdel.f	5.1   11/6/94
c***purpose            
c***description
c   module to assemble the two-dimensional integrals in xyz into
c   primitive one-electron integrals in prmint.
c
c***references
c
c***routines called
c
c***end prologue       fmdel.f
      implicit none
c     --- input variables -----
      integer npint,lenblk
      integer imax,jmax,mini,minj,maxi,maxj,coord
c     --- input arrays (unmodified) ---
      integer nx(*),ny(*),nz(*)
      real*8 xyz(npint,0:imax,0:jmax,3,2)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 prmint(npint,lenblk)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer intgrl,i,ix,iy,iz,j,jx,jy,jz
c
c     --- assemble the one-dimensional integrals in xyz into 
c         three-dimensional integrals in prmint.
      intgrl=0
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
            if(coord.eq.1) then
c              --- (d/dx)*y*z.
               call vmul(prmint(1,intgrl),xyz(1,ix,jx,1,2),
     $                                    xyz(1,iy,jy,2,1),npint)
               call vmul(prmint(1,intgrl),prmint(1,intgrl),
     $                                    xyz(1,iz,jz,3,1),npint)
            else if(coord.eq.2) then
c              --- x*(d/dy)*z.
               call vmul(prmint(1,intgrl),xyz(1,ix,jx,1,1),
     $                                    xyz(1,iy,jy,2,2),npint)
               call vmul(prmint(1,intgrl),prmint(1,intgrl),
     $                                    xyz(1,iz,jz,3,1),npint)
            else if(coord.eq.3) then
c              --- x*y*(d/dz).
               call vmul(prmint(1,intgrl),xyz(1,ix,jx,1,1),
     $                                    xyz(1,iy,jy,2,1),npint)
               call vmul(prmint(1,intgrl),prmint(1,intgrl),
     $                                    xyz(1,iz,jz,3,2),npint)
            endif
    1    continue
    2 continue
c
c
      return
      end
