*deck %W%  %G%
      subroutine prtint(t1,norbs,i,j)
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
      implicit real*8 (a-h,o-z)
      dimension t1(norbs,2)
c
      common /io/ inp,iout
c
c
      write(iout,1) i,j
 1    format(//,' tri. print  row indices ',2(1x,i5))
      write(iout,10)((t1(kk,ll),kk=1,ll),ll=1,i)
 10   format(6(1x,f12.8))
c
      return
      end
