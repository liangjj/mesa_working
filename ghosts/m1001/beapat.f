*deck %W%  %G%
      subroutine beapat(tr,sq,n)
      implicit integer(a-z)
      real*8 tr(n,n),sq(n,n)
c
c
      do 1 i=1,n
         do 2 j=1,n
            tr(i,j)=sq(i,j)+sq(j,i)
    2    continue
    1 continue
c
      return
      end
