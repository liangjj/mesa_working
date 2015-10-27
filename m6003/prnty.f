*deck @(#)prnty.f	1.1 9/8/91
      subroutine prnty(plm,x,n,mu,lmax)
      implicit real *8 (a-h,o-z)
      dimension x(n), plm(n,0:lmax)
      common /io/ inp, iout         
      write (iout,100) mu
      do 200 i=1,n
         write (iout,300) x(i)
         write (iout,400) (plm(i,j),j=mu,lmax)
  200 continue
  100 format(/,5x,'plm for mu',1x,i2)
  300 format(/,5x,'arg',1x,e15.8)
  400 format( (/,5x,5e15.8) )
      return
      end     
