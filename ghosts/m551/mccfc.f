*deck %W%  %G%
      subroutine mccfc(fock,c,fmo,tv,nob1,nob2,nbf1,nbf2)
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
cc
cmp   extended dummy fock,c,tv,fmo
cc
      dimension fock(nbf1,2),c(nbf2,2),tv(2),fmo(nob2,2)
      common / number / zero,pt5,one,two,four,eight
c
c   this program is used to perform a two index
c   transformation on array  fock
c
c
      do 30 i=1,nob1
         do 20 j=1,nbf1
            xx=zero
            do 10 k=1,nbf2
               xx=xx+fock(j,k)*c(k,i)
 10         continue
            tv(j)=xx
 20      continue
         do 40 j=1,nob2
            xx=zero
            do 50 k=1,nbf1
               xx=xx+tv(k)*c(k,j)
 50         continue
            fmo(j,i)=xx
 40      continue
 30   continue
c
      return
      end
