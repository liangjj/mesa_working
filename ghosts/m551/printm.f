*deck %W%  %G%
      subroutine printm(x,len,it)
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
      common /io/ inp,iout
c
cc
cmp   extended dummy x
cc
      dimension x(2)
      if(it.gt.1) go to 15
      write(iout,100)
      ix=0
      do 10 i=1,len
         write(iout,11) i
         write(iout,12) (x(ix+m),m=1,i)
         ix=ix+i
 10   continue
 11   format(/,'  row  ',i6)
 12   format(4(1x,d18.10))
      return
 15   continue
      write(iout,200)
      ix=0
      do 20 i=1,len
         write(iout,11) i
         write(iout,12)(x(ix+m),m=1,len)
         ix=ix+len
 20   continue
c
 100  format(//,'  triangular matrix  ',/)
 200  format(//,'  square     matrix  ',/)
c
      return
      end
