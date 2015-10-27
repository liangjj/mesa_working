*deck %W%  %G%
      subroutine mciiki(dj,dkii,dkki,ni,nk,nki,nii,niis,den,ipair)
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
cmp   extended dummy dj,dkii,dkki,ipair,den
cc
      common /number/ zero,pt5,one,two,four,eight
      dimension dj(nii,2),dkii(nii,2),dkki(nki,2),ipair(2)
      dimension den(2)
c
c----------------------------c
c     mk > mi
c----------------------------c
c
      iiki = 0
      ki=0
c--------------------------------c
c  copy for coulomb matrix
c--------------------------------c
      do 60 i = 1, ni
         do 50 k = 1, nk
            ki=ki+1
            kit=ki+nki
            do 40 ii = 1, nii
               iiki = iiki + 1
               dj(ii,ki) = den(iiki)
               dkii(ii,ki)=zero
               dkii(ii,kit)=zero
 40         continue
            do 45 ii=1,niis
               dkki(ki,ii)=zero
 45         continue
 50      continue
 60   continue
c
c-------------------------------------c
c   build the exchange matrices
c-------------------------------------c
c
      ix=0
      ix1=0
      ii1=0
      do 100 i1=1,ni
         ip1=ipair(i1)
         ix2=0
         ii2=0
         do 90  i2=1,i1
            ip2=ipair(i2)
            ix=ix+1
            ki3=0
            if(i1.eq.i2) go to 85
            do 80 i3=1,ni
               ip3=ipair(i3)
               if(i2.gt.i3) go to 52
               ip32=ip3+i2
               iix1=ix1
               if(i1.gt.i3) go to 51
               ip31=ip3+i1
               iix2=ix2
               go to 53
 51            continue
               ip31=ip1+i3
               iix2=ix2+nki
               go to 53
 52            continue
               ip32=ip2+i3
               ip31=ip1+i3
               iix1=ix1+nki
               iix2=ix2+nki
 53            continue
               do 70 k=1,nk
                  ki3=ki3+1
                  dkki(ix1+k,ii2+i3)=dj(ix,ki3)
                  dkki(ix2+k,ii1+i3)=dj(ix,ki3)
                  dkii(ip31,iix2+k)=dj(ix,ki3)
                  dkii(ip32,iix1+k)=dj(ix,ki3)
 70            continue
 80         continue
c
            go to 185
c
 85         continue
c
            do 180 i3=1,ni
               ip3=ipair(i3)
               if(i1.gt.i3) go to 151
               ip31=ip3+i1
               iix2=ix2
               go to 152
 151           continue
               ip31=ip1+i3
               iix2=ix2+nki
 152           continue
               do 170 k=1,nk
                  ki3=ki3+1
                  dkki(ix1+k,ii2+i3)=two*dj(ix,ki3)
                  dkii(ip31,iix2+k)=two*dj(ix,ki3)
 170           continue
 180        continue
c
 185        continue
c
            ix2=ix2+nk
            ii2=ii2+ni
 90      continue
         ix1=ix1+nk
         ii1=ii1+ni
 100  continue
c
c
c     reorder dkii  for yoshimine convention   ??!!?!"*&
c
c     if(ni.eq.1) return
c
c     ix=1
c     do 995 i=2,ni
c     i1=i-1
c     do 994 j=1,i1
c     ix=ix+1
c     do 990 kj=1,nki
c     xx=dkii(ix,nki+kj)
c     dkii(ix,nki+kj)=dkii(ix,kj)
c     dkii(ix,kj)=xx
c 990 continue
c 994 continue
c     ix=ix+1
c 995 continue
c
      return
      end
