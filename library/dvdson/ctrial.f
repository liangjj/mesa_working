*deck ctrial.f
      subroutine ctrial(mat,buf,ibuf,d,e,e0,v,u0,nel,ntrial,n)
      implicit integer (a-z)
      complex*16 mat, buf, d, e, e0, v, u0
      dimension mat(n,n), buf(*), ibuf(2,*)
      dimension d(*), e(*), e0(*), v(n,*), u0(n,*)
      call czero(mat,n*n)
      do 10 ii=1,nel
         i=ibuf(1,ii)
         j=ibuf(2,ii)
         mat(i,j)=buf(ii)
 10   continue
      do 20 i=1,n
         mat(i,i)=d(i)
 20   continue   
      call rdiag(buf,buf,e,e,e,e,v,v,v,v,v,d,dum,d,n,n,
     1           'real-symmetric')
      return
      end       
