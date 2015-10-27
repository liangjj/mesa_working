*deck %W%  %G%
      subroutine fixbin(xbin,bin,nocc,mij)
C
C***Begin prologue     fixbin
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
C***End prologue       fixbin
C
      implicit real*8 (a-h,o-z)
      dimension bin(2),xbin(2)
c
      ntr=(nocc*(nocc+1))/2
      noc2=nocc*nocc
c
      ix=1
      jx=1
      do 10 i=1,mij
         call totr(xbin(ix),bin(jx),nocc,ntr)
         ix=ix+ntr
         jx=jx+noc2
 10   continue
c
      return
      end
