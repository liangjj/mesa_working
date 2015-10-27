      subroutine mccfc(fock,c,fmo,tv,nob1,nob2,nbf1,nbf2)
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
c     implicit real*8(a-h,o-p,r-z),             integer*2(q)
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
