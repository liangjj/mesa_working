*deck %W%  %G%
      subroutine prtvec(vec,nbf)
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
      real*8 vec(nbf,2)
c
      common /io/ inp,iout
c
c
      write(iout,2)
 2    format(//,'  vector output ')
      do 20 i=1,nbf,6
         lim=min(i+6-1,nbf)
         write(iout,1)(k,k=i,lim)
 1       format(/,6(4x,i3,3x))
cmp   write(iout,3)(k,k=i,lim)
         write(iout,3)
 3       format(/,6(1x,'--------',1x))
         do 10 j=1,nbf
            write(iout,50)(vec(j,m),m=i,lim)
 10      continue
 20   continue
c
 50   format(6(1x,f9.6))
c
      return
      end
