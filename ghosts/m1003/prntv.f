*deck %W%  %G%
      subroutine prntv(v,n)
      implicit real*8(a-h,o-z)
c
       dimension  v(n)
        write(6,1)
  1     format(/,' vector print')
       write(6,2)(v(i),i=1,n)
  2    format(5(1x,f12.8))
c
       return
       end
