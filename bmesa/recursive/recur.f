*deck recur.f
      program recur
c
      implicit integer (a-z)
      real*8 val
      common/io/inp, iout      
      call drum
      n=3
      val=1.d0
      call tstmem(val,n)
      call chainx(0)               
      stop
      end
