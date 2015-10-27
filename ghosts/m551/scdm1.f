*deck %W%  %G%
      subroutine scdm1(dm,nao,nij)
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
      dimension dm(nij)
c
      ix=0
      do 10 i=1,nao
         ix=ix+i
         dm(ix)=dm(ix)*.5d+00
 10   continue
c
      do 20 i=1,nij
         dm(i)=dm(i)*2.d+00
 20   continue
c
      return
      end
