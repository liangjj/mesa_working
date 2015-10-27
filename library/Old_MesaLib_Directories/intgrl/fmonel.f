*deck @(#)fmonel.f	5.1  11/6/94
      subroutine fmonel(prmint,xyz,npint,lenblk,imax,jmax,mini,maxi,
     $                  minj,maxj,nx,ny,nz)
c***begin prologue     fmonel
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)  .
c***keywords           one-electron, integrals
c***author             saxe, paul (lanl).
c***source
c***purpose            assembles the two-dimensional integrals into primitive
c                      one-electron integrals.
c***description
c                      call fmonel(prmint,xyz,npint,lenblk,imax,jmax,mini,
c                                  maxi,minj,maxj,nx,ny,nz)
c
c                        prmint   sink for primitive integrals.
c                        xyz      source of two-dimensional integrals.
c                        npint    number of primitive integrals to form.
c                        lenblk
c
c***references
c***routines called    vmul(math)
c***end prologue       fmonel
      implicit integer (a-z)
c
      real*8 prmint(npint,lenblk),xyz(npint,0:imax,0:jmax,3)
      integer nx(*),ny(*),nz(*)
c
c     ----- start timing -----
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
            call vmul(prmint(1,intgrl),xyz(1,ix,jx,1),xyz(1,iy,jy,2),
     #                                                         npint)
            call vmul(prmint(1,intgrl),prmint(1,intgrl),xyz(1,iz,jz,3),
     #                                                         npint)
    1    continue
    2 continue
c
c     ----- stop timing -----
c
c
c
      return
      end
