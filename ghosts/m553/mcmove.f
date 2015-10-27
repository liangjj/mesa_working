      subroutine mcmove(c,ct,nob,nbf)
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
