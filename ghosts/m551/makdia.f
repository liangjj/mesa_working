*deck %W%  %G%
      subroutine makdia(buf,diag,n)
c
c***begin prologue     makdia
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
c***end prologue       makdia
c
      implicit real*8(a-h,o-z)
c
      dimension buf(n,n),diag(2)
c
      ix=0
      do 10 i=1,n
         diag(i)=buf(i,i)
 10   continue
c
      return
      end
