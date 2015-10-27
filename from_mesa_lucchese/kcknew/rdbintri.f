      subroutine rdbintri(amat,nrow,iunit)
       implicit real*8 (a-h,o-z)
      dimension amat(nrow*(nrow+1)/2)
      nn=nrow*(nrow+1)/2
      read(iunit) (amat(i),i=1,nn)
      return
      end
