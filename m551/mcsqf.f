*deck @(#)mcsqf.f	1.1  11/30/90
      subroutine mcsqf(fab,f,nbf)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcsqf.f	1.1   11/30/90
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
cmp   extended dummy fab,f
cc
      dimension fab(2),f(nbf,nbf)
c
      ix=0
      do 10 i=1,nbf
         do 20 j=1,i
            ix=ix+1
            f(i,j)=fab(ix)
            f(j,i)=fab(ix)
 20      continue
 10   continue
c
      return
      end
