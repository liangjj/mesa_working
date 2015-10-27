      subroutine vamul(a,b,c,n)
      implicit real*8(a-h,o-z)
c
      dimension a(n),b(n),c(n)
c
      do 10 i=1,n
       a(i)=a(i)+b(i)*c(i)
  10  continue
c
      return
      end
