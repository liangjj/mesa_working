*deck mul3x3
      subroutine mul3x3(a,b,c)
      implicit real*8(a-h,o-z)
c
c     a little routine to do a 3x3 matrix multtiplication.
c     returns in "c" the product "ab".
c
      dimension a(3,3),b(3,3),c(3,3)
      data zero/0.d0/
c
c     call rtrace(6hmul3x3,1)
      do 100 i=1,3
         do 100 j=1,3
            c(i,j)=zero
            do 100 k=1,3
 100           c(i,j)=c(i,j)+a(i,k)*b(k,j)
               return
               end
