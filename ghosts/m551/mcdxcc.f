*deck %W%  %G%
      subroutine mcdxcc(dab,hess,xjk,nab,mrs,itran,thrsh)
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
cc
cmp   extended dummy dab,hess,xjk
cc
      dimension dab(2),hess(2),xjk(2)
c
      if(itran.ne.0)go to 200
c
      nx=0
      do 120 i=1,nab
         den=dab(i)
         if(abs(den).lt.thrsh)go to 115
         do 110 j=1,mrs
            hess(nx+j)=hess(nx+j)+den*xjk(j)
 110     continue
 115     continue
         nx=nx+mrs
 120  continue
c
      return
c
 200  continue
c
      do 220 j=1,mrs
         xint=xjk(j)
         if(abs(xint).lt.thrsh)go to 220
         nx=j
         do 210 i=1,nab
            hess(nx)=hess(nx)+xint*dab(i)
            nx=nx+mrs
 210     continue
 220  continue
c
      return
      end
