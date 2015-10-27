*deck %W%  %G%
      subroutine movham(ham,temp,n)
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
cmp   extended dummy ham,temp
      dimension ham(2),temp(2)
ccccc
c     copy lower tria. matrix
cccc
      ii=0
      do 20 i=1,n
         do 10 j=1,i
            ii=ii+1
            temp(ii)=ham(ii)
 10      continue
 20   continue
c
      return
      end
