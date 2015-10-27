      subroutine wrbinsqr(amat,nrow,ncol,iunit)
       implicit real*8 (a-h,o-z)
      dimension amat(nrow,ncol)
      do 1 i=1,ncol
  1   write(iunit) (amat(j,i),j=1,ncol)
      return
      end
