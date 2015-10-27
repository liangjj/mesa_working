*deck rfill.f
      subroutine rfill(d,e,e0,v,u0,n)
      implicit integer (a-z)
      real*8 d, e, e0, v, u0
      dimension d(*), e(*), e0(*), v(n,*), u0(n,*)
      do 10 i=1,n
         e(i)=d(i)
         e0(i)=d(i)
         v(i,i)=1.d0
         u0(i,i)=1.d0
 10      continue
      return
      end       
