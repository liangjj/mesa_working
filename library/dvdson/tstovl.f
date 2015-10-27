*deck tstovl.f
      subroutine tstovl(vec,cvec,n,nvc,mattyp)
      implicit integer (a-z)
      character*(*) mattyp
      real*8 vec, ovl, sdot
      complex*16 cvec, covl, cdotc
      dimension vec(n,nvc), cvec(n,nvc)
      common/io/inp, iout
      if(mattyp.eq.'complex') then
         do 10 i=1,nvc
            do 20 j=1,i
               covl = cdotc(n,cvec(1,i),1,cvec(1,j),1)
               write(iout,1) i, j, covl
   20       continue
   10    continue
      else
         do 30 i=1,nvc
            do 40 j=1,i
               ovl = sdot(n,vec(1,i),1,vec(1,j),1)
               write(iout,2) i, j, ovl
   40       continue
   30    continue
      endif          
      return
    1 format(1x,'i = ',i3,1x,'j = ',i3,1x,'overlap = ',e15.8,1x,e15.8)
    2 format(1x,'i = ',i3,1x,'j = ',i3,1x,'overlap = ',e15.8)
      end       
