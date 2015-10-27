*deck %W%  %G%
      subroutine tranaa(xm,c,s,r,nob,nbf)
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
c
c      extended dummy c,xm,s,r
c
      dimension c(nbf,nob),xm(nbf,nbf),s(nbf,nob),r(2)
c
      do 30 i=1,nbf
         do 20 j=1,nob
            xx=0.d0
            do 15 k=1,nbf
               xx=xx+xm(i,k)*c(k,j)
 15         continue
            s(i,j)=xx
 20      continue
 30   continue
c
c
      ix=0
      do 100 i=1,nob
         do 90  j=1,i
            xx=0.d0
            do 80 k=1,nbf
               xx=xx+c(k,i)*s(k,j)
 80         continue
            ix=ix+1
            r(ix)=xx
 90      continue
 100  continue
c
      return
      end
