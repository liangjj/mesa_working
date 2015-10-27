*deck @(#)fmt.f	5.1  11/6/94
      subroutine fmt(prmint,xyz,npint,lenblk,imax,jmax,t1,mini,maxi,
     #               minj,maxj,nx,ny,nz,nderiv,nint)
c
c***module to assemble the two-dimensional integrals in xyz into
c   primitive kinetic-energy integrals in prmint.
c
c paul saxe                      23 july 1984                lanl
c
      implicit integer (a-z)
c
      real*8 prmint(npint,lenblk,nint)
      real*8 xyz(npint,0:imax,0:jmax,3,0:nderiv,2)
      real*8 t1(npint)
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
            do 10 k=1,npint
               prmint(k,intgrl,1)=xyz(k,ix,jx,1,0,2)*
     #             (xyz(k,iy,jy,2,0,1)*xyz(k,iz,jz,3,0,1))+
     #                     xyz(k,ix,jx,1,0,1)*
     #             (xyz(k,iy,jy,2,0,2)*xyz(k,iz,jz,3,0,1)+
     #              xyz(k,iy,jy,2,0,1)*xyz(k,iz,jz,3,0,2))
   10       continue
            if (nderiv.ge.1) then
               do 11 k=1,npint
                  prmint(k,intgrl,2)=xyz(k,ix,jx,1,1,2)*
     #                (xyz(k,iy,jy,2,0,1)*xyz(k,iz,jz,3,0,1))+
     #                        xyz(k,ix,jx,1,1,1)*
     #                (xyz(k,iy,jy,2,0,2)*xyz(k,iz,jz,3,0,1)+
     #                 xyz(k,iy,jy,2,0,1)*xyz(k,iz,jz,3,0,2))
                  prmint(k,intgrl,3)=xyz(k,ix,jx,1,0,2)*
     #                (xyz(k,iy,jy,2,1,1)*xyz(k,iz,jz,3,0,1))+
     #                        xyz(k,ix,jx,1,0,1)*
     #                (xyz(k,iy,jy,2,1,2)*xyz(k,iz,jz,3,0,1)+
     #                 xyz(k,iy,jy,2,1,1)*xyz(k,iz,jz,3,0,2))
                  prmint(k,intgrl,4)=xyz(k,ix,jx,1,0,2)*
     #                (xyz(k,iy,jy,2,0,1)*xyz(k,iz,jz,3,1,1))+
     #                        xyz(k,ix,jx,1,0,1)*
     #                (xyz(k,iy,jy,2,0,2)*xyz(k,iz,jz,3,1,1)+
     #                 xyz(k,iy,jy,2,0,1)*xyz(k,iz,jz,3,1,2))
   11          continue
            end if
            if (nderiv.ge.2) then
               do 12 k=1,npint
                  prmint(k,intgrl,5)=xyz(k,ix,jx,1,2,2)*
     #                (xyz(k,iy,jy,2,0,1)*xyz(k,iz,jz,3,0,1))+
     #                        xyz(k,ix,jx,1,2,1)*
     #                (xyz(k,iy,jy,2,0,2)*xyz(k,iz,jz,3,0,1)+
     #                 xyz(k,iy,jy,2,0,1)*xyz(k,iz,jz,3,0,2))
                  prmint(k,intgrl,6)=xyz(k,ix,jx,1,1,2)*
     #                (xyz(k,iy,jy,2,1,1)*xyz(k,iz,jz,3,0,1))+
     #                        xyz(k,ix,jx,1,1,1)*
     #                (xyz(k,iy,jy,2,1,2)*xyz(k,iz,jz,3,0,1)+
     #                 xyz(k,iy,jy,2,1,1)*xyz(k,iz,jz,3,0,2))
                  prmint(k,intgrl,7)=xyz(k,ix,jx,1,0,2)*
     #                (xyz(k,iy,jy,2,2,1)*xyz(k,iz,jz,3,0,1))+
     #                        xyz(k,ix,jx,1,0,1)*
     #                (xyz(k,iy,jy,2,2,2)*xyz(k,iz,jz,3,0,1)+
     #                 xyz(k,iy,jy,2,2,1)*xyz(k,iz,jz,3,0,2))
                  prmint(k,intgrl,8)=xyz(k,ix,jx,1,1,2)*
     #                (xyz(k,iy,jy,2,0,1)*xyz(k,iz,jz,3,1,1))+
     #                        xyz(k,ix,jx,1,1,1)*
     #                (xyz(k,iy,jy,2,0,2)*xyz(k,iz,jz,3,1,1)+
     #                 xyz(k,iy,jy,2,0,1)*xyz(k,iz,jz,3,1,2))
                  prmint(k,intgrl,9)=xyz(k,ix,jx,1,0,2)*
     #                (xyz(k,iy,jy,2,1,1)*xyz(k,iz,jz,3,1,1))+
     #                        xyz(k,ix,jx,1,0,1)*
     #                (xyz(k,iy,jy,2,1,2)*xyz(k,iz,jz,3,1,1)+
     #                 xyz(k,iy,jy,2,1,1)*xyz(k,iz,jz,3,1,2))
                  prmint(k,intgrl,10)=xyz(k,ix,jx,1,0,2)*
     #                (xyz(k,iy,jy,2,0,1)*xyz(k,iz,jz,3,2,1))+
     #                        xyz(k,ix,jx,1,0,1)*
     #                (xyz(k,iy,jy,2,0,2)*xyz(k,iz,jz,3,2,1)+
     #                 xyz(k,iy,jy,2,0,1)*xyz(k,iz,jz,3,2,2))
   12          continue
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
