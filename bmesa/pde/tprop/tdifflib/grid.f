      subroutine grid(x,x0,x1,h,nstp)
      implicit integer(a-z)
      real*8 x, x0, x1, val, h
      character*80 title
      dimension x(*)
      common/io/inp, iout
      x(1)=x0
      val=x(1)
      nstp=1
      do while (val.lt.x1)
         val=val+h
         nstp=nstp+1
         x(nstp)=val
      enddo
      if(x(nstp).gt.x1) then
         nstp=nstp-1
      endif
c      title='grid'
c      call prntrm(title,x,nstp,1,nstp,1,iout)
      return
      end
