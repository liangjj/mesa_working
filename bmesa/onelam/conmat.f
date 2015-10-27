      subroutine conmat(s1,s2,v,lamat,wt,n)
      implicit integer (a-z)
      real *8 s1, s2, v, lamat, wt
      dimension s1(n), s2(n), v(n), lamat(n,n), wt(n)
      call rzero(lamat,n*n)
      do 10 i=1,n
         lamat(i,i)=1.e+00
   10 continue
      do 20 i=1,n
         do 30 j=1,i
            lamat(i,j)=lamat(i,j)-s2(i)*s1(j)*v(j)*wt(j)
   30    continue
   20 continue
      do 40 i=1,n
         do 50 j=i+1,n
            lamat(i,j)=lamat(i,j)-s2(j)*s1(i)*v(j)*wt(j)
   50    continue
   40 continue
      return
      end
