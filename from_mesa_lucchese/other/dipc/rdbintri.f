      subroutine rdbintri(amat,nrow,iunit)
      real*8 amat(nrow*(nrow+1)/2)
      nn=nrow*(nrow+1)/2
      read(iunit) (amat(i),i=1,nn)
      return
      end
