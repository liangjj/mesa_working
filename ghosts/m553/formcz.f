      subroutine formcz(cz,nbf)
c
c***begin prologue     formcz
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
c***end prologue       formcz
c
      implicit real*8 (a-h,o-z)
c
      dimension cz(nbf,2)
c
      common / number / zero,pt5,one,two,four,eight
c
      do 30 i=1,nbf
         do 20 j=1,nbf
            cz(i,j)=zero
 20      continue
         cz(i,i)=one
 30   continue
c
      return
      end
