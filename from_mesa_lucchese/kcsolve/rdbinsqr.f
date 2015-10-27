      subroutine rdbinsqr(amat,nrow,ncol,iunit)
      real*8 amat(nrow,ncol)
c      write (6,*)'rdbinsqr nrow', nrow, '  ncol', ncol
c      write(6,*) iunit, 'iunit'
      do 1 i=1,ncol
c         write (6,*) 'icol', i
  1   read(iunit) (amat(j,i),j=1,ncol)
      return
      end
