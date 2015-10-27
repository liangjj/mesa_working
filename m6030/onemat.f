      subroutine onemat(smat,ncon)
      implicit integer (a-z)
      real *8 smat
      dimension smat(ncon,ncon)
      call rzero(smat,ncon*ncon)
      do 10 i=1,ncon
         smat(i,i)=1.e+00
   10 continue
      return
      end
