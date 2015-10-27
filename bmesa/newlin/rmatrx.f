      subroutine rmatrx (rmat,rhs,ncst,nsts,nsol,ntchn,nptmx,msize,ien
     1 ,ops)
      implicit integer(a-z)
      common /io/ inp, iout
      character *(*) ops
      character *3 itoc
      real *8rmat, rhs
      dimension rhs(msize,nsol), rmat(ntchn,ntchn), ncst(nsts)
*
      do 10 i=1,ntchn
      icnt=0
      jj=0
      do 10 j=1,nsts
      nn=ncst(j)
      do 10 k=1,nn
      icnt=icnt+1
      jj=jj+nptmx
      rmat(i,icnt)=-rhs(jj,i)
   10 continue
      call iosys ('write real "r matrix'//itoc(ien)//'" to rwf',ntchn
     1 *ntchn,rmat,0,0)
      write (iout,20)
      call matprt (rmat,ntchn,ntchn,ntchn,ntchn,0,0,0,0,0,0,0)
*
      return
c
   20 format (//,20x,'r - matrix',/)
      end
