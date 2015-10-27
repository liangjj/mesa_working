*deck @(#)mczero.f	1.1  11/30/90
      subroutine mczero(f,nbf,nob)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mczero.f	1.1   11/30/90
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
cmp   extended dummy f
cc
      common / number / zero,pt5,one,two,four,eight
      dimension f(2)
c
c
      ix=0
      do 11 j=1,nob
         do 10 i=1,nbf
            ix=ix+1
            f(ix)=zero
 10      continue
 11   continue
c
      return
      end
