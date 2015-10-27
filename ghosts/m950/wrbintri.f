*deck @(#)wrbintri.f	1.1  11/30/90
      subroutine wrbintri(amat,nrow,iunit)
      dimension amat(nrow*(nrow+1)/2)
      nn=nrow*(nrow+1)/2
      write(iunit) (amat(i),i=1,nn)
      return
      end
