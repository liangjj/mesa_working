      subroutine cmwrite(h,n,m,nmax)
      implicit real*8(a-h,o-z)
      complex*16 h(nmax,nmax)
      do 1 i=1,n
 1    write(6,2)i,(h(i,j),j=1,m)
 2    format(i5,8e12.4/(5x,8e12.4))
      return
      end

