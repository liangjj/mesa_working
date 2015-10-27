      function dotpro(n,x,ix,y,iy)
      dimension x(ix,n),y(iy,n)
      dotpro=0.0
      if(n .le. 0 ) go to 20
      do 10 i=1,n
      dotpro=dotpro+x(1,i)*y(1,i)
 10   continue
 20   return
      end
