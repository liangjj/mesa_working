*deck %W%  %G%
      subroutine mcmove(c,ct,nob,nbf)
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
c      extended dummy c,ct
c
      dimension c(nbf,nob),ct(nbf,nob)
c
      do 10 i=1,nob
         do 20 j=1,nbf
            ct(j,i)=c(j,i)
 20      continue
 10   continue
c
      return
      end
