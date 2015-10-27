      subroutine vecprod(b,x,a,z,n)
      dimension x(n),z(n)
      if(n .le. 0 ) go to 20
      do 10 i=1,n
      z(i)=a*x(i)+b
 10   continue
 20   return
      end
