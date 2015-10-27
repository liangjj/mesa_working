      subroutine zimag(array,n,m)
      real *8 r1
      complex *16 array(n,m)
      do 10 i=1,n
         do 20 j=1,m
            r1=real(array(i,j))
            array(i,j)=dcmplx(r1,0.d0)
   20    continue
   10 continue  
      return
      end
