*deck cc2opy
      subroutine cc2opy(a,b,n)
      implicit integer (a-z)
      complex*16 a,b
      dimension a(n), b(n)
      do 10 i=1,n
   10 b(i)=a(i)
      return
      end
