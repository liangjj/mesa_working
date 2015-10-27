*deck tstnrm
      subroutine tstnrm(ply,wt,smat,nq,nprim)
      implicit integer (a-z)
      real *8 ply, smat, wt
      dimension ply(nq,nprim), smat(nprim,nprim), wt(nq)
      character*80 title
      common/io/ inp, iout
c     overlap matrix elements
c
      do 10 i=1,nprim
         do 20 j=1,i
            smat(i,j)=0.d0
            do 30 k=1,nq
               smat(i,j)=smat(i,j)+ply(k,i)*wt(k)*ply(k,j)
   30       continue
            smat(j,i)=smat(i,j)
   20    continue
   10 continue
      title='overlap matrix'
      call prntrm(title,smat,nprim,nprim,nprim,nprim,iout) 
      return
      end
