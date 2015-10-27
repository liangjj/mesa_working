*deck @(#)wrbinsqr.f	1.1  11/30/90
c
      subroutine wrbinsqr(amat,nrow,ncol,iunit)
      dimension amat(nrow,ncol)
      do 1 i=1,ncol
  1   write(iunit) (amat(j,i),j=1,ncol)
      return
      end
