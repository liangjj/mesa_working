*deck initv
      subroutine initv(v,n)
      implicit integer(a-z)
      real*8 v
      dimension v(n)
      call rzero(v,n)
      do 10 i=1,n
         v(i)=1.d0
 10   continue   
      return
      end
