*deck @(#)fmoned.f	5.1  11/6/94
      subroutine fmoned(prmint,xyz,npint,lenblk,imax,jmax,mini,maxi,
     #                  minj,maxj,nx,ny,nz,nderiv,nint)
c
c***module to assemble the two-dimensional integrals in xyz into
c   primitive one-electron integrals in prmint.
c
c paul saxe                      23 july 1984                lanl
c
      implicit integer (a-z)
c
      real*8 prmint(npint,lenblk,nint)
      real*8 xyz(npint,0:imax,0:jmax,3,0:nderiv)
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
            call vmul(prmint(1,intgrl,1),xyz(1,ix,jx,1,0),
     #                             xyz(1,iy,jy,2,0),npint)
            call vmul(prmint(1,intgrl,1),prmint(1,intgrl,1),
     #                             xyz(1,iz,jz,3,0),npint)
            if (nderiv.ge.1) then
               do 55 k=1,npint
                  prmint(k,intgrl,2)=xyz(k,ix,jx,1,1)*
     #                               xyz(k,iy,jy,2,0)*
     #                               xyz(k,iz,jz,3,0)
                  prmint(k,intgrl,3)=xyz(k,ix,jx,1,0)*
     #                               xyz(k,iy,jy,2,1)*
     #                               xyz(k,iz,jz,3,0)
                  prmint(k,intgrl,4)=xyz(k,ix,jx,1,0)*
     #                               xyz(k,iy,jy,2,0)*
     #                               xyz(k,iz,jz,3,1)
   55          continue
            end if
            if (nderiv.ge.2) then
               do 56 k=1,npint
                  prmint(k,intgrl,5)=xyz(k,ix,jx,1,2)*
     #                               xyz(k,iy,jy,2,0)*
     #                               xyz(k,iz,jz,3,0)
                  prmint(k,intgrl,6)=xyz(k,ix,jx,1,1)*
     #                               xyz(k,iy,jy,2,1)*
     #                               xyz(k,iz,jz,3,0)
                  prmint(k,intgrl,7)=xyz(k,ix,jx,1,0)*
     #                               xyz(k,iy,jy,2,2)*
     #                               xyz(k,iz,jz,3,0)
                  prmint(k,intgrl,8)=xyz(k,ix,jx,1,1)*
     #                               xyz(k,iy,jy,2,0)*
     #                               xyz(k,iz,jz,3,1)
                  prmint(k,intgrl,9)=xyz(k,ix,jx,1,0)*
     #                               xyz(k,iy,jy,2,1)*
     #                               xyz(k,iz,jz,3,1)
                  prmint(k,intgrl,10)=xyz(k,ix,jx,1,0)*
     #                               xyz(k,iy,jy,2,0)*
     #                               xyz(k,iz,jz,3,2)
   56          continue
            end if
    1    continue
    2 continue
c
c     ----- stop timing -----
c
c
c
      return
      end
