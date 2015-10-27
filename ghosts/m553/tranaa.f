      subroutine tranaa(xm,c,s,r,nob,nbf)
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
