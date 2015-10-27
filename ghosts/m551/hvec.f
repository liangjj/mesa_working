*deck %W%  %G%
      subroutine hvec(buf,c,b,n,nvec)
c
c***begin prologue     hvec
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
c***end prologue       hvec
c
      implicit real*8(a-h,o-z)
c
      dimension buf(n,n),c(n,nvec),b(n,nvec)
c
      call ebc(b,buf,c,n,n,nvec)
c
      return
      end
