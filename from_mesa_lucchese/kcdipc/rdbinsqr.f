      subroutine rdbinsqr(amat,nrow,ncol,iunit)
      real*8 amat(nrow,ncol)
      do 1 i=1,ncol
  1   read(iunit) (amat(j,i),j=1,ncol)
      return
      end
