      subroutine vscale(vec,fac,n)
      implicit integer (a-z)
      real *8 vec, fac
      dimension vec(n)
      do 10 i=1,n
         vec(i)=vec(i)*fac
   10 continue
      return
      end
