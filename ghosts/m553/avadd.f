      subroutine avadd(a,b,n)
c
c***begin prologue     avadd
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
c***end prologue       avadd
c
      dimension a(n),b(n)
c
      do 1 i=1,n
         a(i)=a(i)+b(i)
 1    continue
c
      return
      end
