*deck %W%  %G%
      subroutine mcstv(result,vector,scalar,length)
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
cc
cmp   extended dummy result,vector
cc
      dimension result(2),vector(2)
c
      do 10 i=1,length
         result(i)=result(i)+scalar*vector(i)
 10   continue
c
      return
      end
