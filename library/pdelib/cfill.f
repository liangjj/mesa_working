*deck cfill.f
      subroutine cfill(d,e,e0,v,u0,n)
      implicit integer (a-z)
      complex*16 d, e, e0, v, u0
      dimension d(*), e(*), e0(*), v(n,*), u0(n,*)
      do 10 i=1,n
         e(i)=d(i)
         e0(i)=d(i)
         v(i,i)=(1.d0,0.d0)
         u0(i,i)=(1.d0,0.d0)
 10      continue
      return
      end       
