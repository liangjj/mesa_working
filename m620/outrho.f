*deck @(#)outrho.f	5.1  11/6/94
      subroutine outrho(amat,nq)
      implicit none
c
      integer nq
      real*8 amat(1)
c
      integer len,jk,i,j,jkend
      integer inp,iout
      common/io/inp,iout
c
  100 format(' correlation matrix, upper triangular form')
  110 format(10f8.3)
  120 format(' ')
c
      write(iout,100)
c
      len=nq
      jk=nq+1
      do 10 j=1,nq
         jkend=jk+len-1
         write(iout,110) (amat(i),i=jk,jkend)
         write(iout,120)
         jk=jk+len
         len=len-1
   10 continue
c
c
      return
      end
