*deck %W%  %G%
      subroutine copy(a,b,n)
      implicit integer (a-z)
      real*8 a,b
      dimension a(n), b(n)
      do 10 i=1,n
         b(i)=a(i)
   10 continue
      return
      end
