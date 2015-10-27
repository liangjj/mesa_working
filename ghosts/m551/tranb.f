*deck %W%  %G%
      subroutine tranb(cv,ct,cp,mobs,nbfs)
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
c      extended dummy cv,ct,cp
c
      dimension cv(2),cp(2),ct(2)
c
      mtot=mobs*nbfs
      do 105 k=1,mtot
         cp(k)=0.d0
 105  continue
c
      kl=0
      koff=0
      do 130 k=1,mobs
         loff=0
         do 125 l=1,mobs
            kl=kl+1
            xx=cv(kl)
            do 120 m=1,nbfs
               cp(koff+m)=cp(koff+m)+xx*ct(loff+m)
 120        continue
            loff=loff+nbfs
 125     continue
         koff=koff+nbfs
 130  continue
c
      return
      end
