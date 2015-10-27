*deck %W%  %G%
      subroutine mciijk(dj,dk,ni,nj,nk,nii,njk,nij,nik,den,ipair)
C
C***Begin prologue
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
C
C***References
C
C***Routines called    (none)
C
C***End prologue
C
      implicit real*8(a-h,o-z)
cc
cmp   extended dummy dj,dk,den,ipair
cc
      dimension dj(nii,njk),dk(nij,nik)
      dimension den(2),ipair(2)
      common /number/ zero,pt5,one,two,four,eight
c
c---------------------------------c
c     mj > mi > mk
c---------------------------------c
c
      iijk=0
      jk=0
      kx=0
      do 40 k=1,nk
         do 30 j=1,nj
            jk=jk+1
            ix=0
            ji=j
            do 20 i = 1, ni
               jii=j
               do 10 ii = 1, i
                  ix=ix+1
                  iijk=iijk+1
                  dj(ix,jk) = den(iijk)
                  dk(ji,kx+ii) = den(iijk)
                  dk(jii,kx+i) = den(iijk)
                  jii=jii+nj
 10            continue
               dk(ji,kx+i)= two*dk(ji,kx+i)
               ji=ji+nj
 20         continue
 30      continue
         kx=kx+ni
 40   continue
c
      return
      end
