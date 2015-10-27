*deck %W%  %G%
      subroutine mcdmsq(den,dab,nob)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
      implicit real*8(a-h,o-z)
cc
cmp   extended dummy den,dab
cc
      dimension den(2),dab(nob,2)
c
      common /io/ inp,iout
c
c
c-------------------------------------------------------c
c      this routine squares the one electron density    c
c       and doubles the diagonal elements               c
c-------------------------------------------------------c
c
      common / number / zero,pt5,one,two,four,eight
c
c     write(iout,7001) nob
c7001 format('  nob  mcdmsq ',i8)
      ix=0
      do 20 i=1,nob
         do 10 j=1,i
            ix=ix+1
            dab(i,j)=den(ix)
            dab(j,i)=den(ix)
 10      continue
         dab(i,i)=dab(i,i)+dab(i,i)
 20   continue
cccccc
c     write(iout,21)
c  21 format('  dab  in mcdmsq ')
c     call printr(dab,nob,nob)
cccccc
      return
      end
