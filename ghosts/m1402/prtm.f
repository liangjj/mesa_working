*deck %W%  %G%
      subroutine prtm(a,n)
      implicit real*8(a-h,o-z)
      dimension a(*)
c
      common /io/ inp,iout
c
   1  format(5(1x,f12.8))
   2  format(/,' row = ',i5)
      ix=0
      do 10 i=1,n
      write(iout,2) i
      write(iout,1)(a(ix+j),j=1,i)
      ix=ix+i
  10  continue
c
      return
      end
