*deck %W%  %G%
      subroutine printi(ipt,len,ii)
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
      implicit real*8 (a-h,o-z)
      dimension ipt(*)
c
      common /io/ inp,iout
c
c
      ix=0
      do 1 i=1,len
         write(iout,2)(ipt(ix+j),j=1,i)
         ix=ix+i
 1    continue
 2    format(5(2x,i8))
c
      return
      end
