*deck @(#)wvec.f	5.1  11/6/94
      subroutine wvec(a,eig,nrows,ncols,rowlab,collab)
c***begin prologue     wvec.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             unknown
c***source             @(#)wvec.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       wvec.f
      implicit integer(a-z)
      character*(*) rowlab(nrows), collab(ncols)
      real*8 a(nrows,ncols), eig(ncols)
c
c     a eigenvector output routine which deals with row and/or column
c     headings.
c
c        a      ... the matrix to be printed, real(nrows,ncols).
c        nrows  ... the number of rows of a to print.
c        ncols  ... the number of columns of a to print.
c        lrow   ... determines kind of row labels desired.
c                   0 ... numbers only.
c                   1 ... label only.
c                   2 ... number and label.
c        lcol   ... determines kind of column labels desired.
c                   0 ... numbers only.
c                   1 ... label only.
c                   2 ... number and label.
c        rowlab ... the row labels. character*(*) (nrows).
c        collab ... the column labels. character*(*) (ncols).
c        isym   ... the type of matrix.
c                   0 ... a full matrix.
c                   1 ... only print lower triangle.
c                  -n ... only one column of the matrix is printed.
c                   the abs(isym)-th column will appear as a row in the
c                   output. note that the column labels in the calling list
c                   refer to the columns of the output, thus in this case
c                   to the row of the matrix.
c        eig    ... the eigenvalues, real(ncols).
c
c
c
c
      call matprt(a,nrows,ncols,nrows,ncols,1,0,rowlab,collab,0,eig,
     $           .true.)
c
c
      return
      end