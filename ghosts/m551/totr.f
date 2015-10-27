*deck %W%  %G%
      subroutine totr(triang,square,num,nnp)
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
c
      implicit integer (a-z)
c
      real*8 triang(nnp),square(num,num)
c
      ij=0
      do 2 i=1,num
         do 1 j=1,i
            ij=ij+1
            triang(ij)=square(j,i)
    1    continue
    2 continue
c
      return
      end
