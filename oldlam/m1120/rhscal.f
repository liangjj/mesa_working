      subroutine rhscal (gr,glast,rhs,nptmx)
      implicit integer(a-z)
      real *8 gr, glast, rhs
      dimension gr(nptmx), rhs(nptmx)
      common /io/ inp, iout
      call rzero(rhs,nptmx)
      do 10 ip=1,nptmx
         rhs(ip)=gr(ip)*glast
   10 continue
      return
c
      end
