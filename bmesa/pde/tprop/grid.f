      subroutine grid(x,x0,h,nstp)
      implicit integer(a-z)
      real*8 x, x0, h
      character*80 title
      dimension x(*)
      common/io/inp, iout
      x(1)=x0
      do 10 i=2,nstp
         x(i)=x(i-1)+h
 10   continue
      title='grid'
      call prntrm(title,x,nstp,1,nstp,1,iout)
      return
      end
