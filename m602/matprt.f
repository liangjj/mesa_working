*deck @(#)matprt.f	1.2  7/30/91
      subroutine matprt(a,md,nd,nrows,ncols,lrow,lcol,rowlab,
     1                  collab,isym,eig,ifeig)
      implicit integer(a-z)
      character*(*) rowlab(nrows), collab(ncols)
      character line*80,tmpcol*8,tmprow*16
      logical ifeig
      real*8 a(md,nd), eig(ncols)
c
      common/io/inp,iout
      parameter (coldat=5)
c
 9000 format(26x,10(4x,i6))
 9010 format(10x,'eigenvalues -- ',1x,10f10.5)
 9020 format(1x,a80)
 9030 format(1x,i9,11x,5f10.5)
 9040 format(1x,i4,1x,10x,a8,2x,5f10.5)
 9050 format(1x,i4,1x,2x,a16,2x,5f10.5)
 9060 format(1x,i4,1x,4x,i3,2x,a8,3x,5f10.5)
 9070 format(1x,i4,1x,i3,1x,a16,5f10.5)
c
c     a matrix output routine to deal with possible row and/or column
c     headings.
c
c        a      ... the matrix to be printed, real(md,nd).
c                   it is filled to (nrows,ncols).
c        md     ... leading dimension of a.
c        nd     ... column dimension of a.
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
c        ifeig  ... whether to print eigenvalues, a logical variable.
c
c
c     for row labels, the routine figures out how wide the label field
c     must be (tabs to 8 or 16 characters).  for column labels, the
c     maximum number of characters is 8.
c
c
c     initialization.
      rowmin=1
      rowmax=nrows
      colmin=0
      colmax=0
c
c     scan the row labels, find the longest.
      maxlen=0
      rowwid=8
      if(lrow.ne.0) then
         do 10 i=1,nrows
            lenr=cskipb(rowlab(i),' ')+1
            maxlen=max(maxlen,lenr)
   10    continue
      endif
      if(maxlen.gt.8) rowwid=16
c
c     top of loop over "passes"; it may take more than one
c     line to get an entire row printed.
      npass=(ncols-1)/coldat+1
      if(isym.lt.0) npass=(nrows-1)/coldat+1
c
      do 100 ipass=1,npass
c        get new column limits.
         colmin=colmax+1
         colmax=colmax+coldat
         colmax=min(ncols,colmax)
         if(isym.lt.0) then
            rowmax=-isym
         else if(isym.eq.1) then
            rowmin=colmin
         endif
c
c        print heading, integers first, then column labels.
         if(lcol.ne.1) write(iout,9000) (i,i=colmin,colmax)
c
         if(lcol.ne.0) then
            line=' '
            ltab=22
            do 20 j=colmin,colmax
               lenc=min(len(collab(j)),8)
               call crjust(collab(j)(1:lenc),tmpcol(1:lenc))
               line(ltab+1:)=tmpcol(1:lenc)
               ltab=ltab+10
   20       continue
            write(iout,9020) line
         endif
c
c        write eigenvalues if requested.
         if(ifeig) write(iout,9010) (eig(i),i=colmin,colmax)
c
c        loop over rows of output.
         do 90 row=rowmin,rowmax
            if(lrow.ne.0) then
               line=' '
               lenr=min(rowwid,len(rowlab(row)))
               call crjust(rowlab(row)(1:lenr),tmprow(1:lenr))
               line(1:lenr)=tmprow(1:lenr)
            endif
c           deal with the possibility of a lower triangle.
            collim=colmax
            if(isym.gt.0) collim=min(row,colmax)
c
c           write the row.
            if(lrow.eq.0) then
               write(iout,9030) row,(a(row,j),j=colmin,collim)
            else if(lrow.eq.1) then
               if(rowwid.eq.8)  write(iout,9040) row,line,
     1                               (a(row,j),j=colmin,collim)
               if(rowwid.eq.16) write(iout,9050) row,line,
     1                               (a(row,j),j=colmin,collim)
            else if(lrow.eq.2) then
               if(rowwid.eq.8)  write(iout,9060) row,line,
     1                               (a(row,j),j=colmin,collim)
               if(rowwid.eq.16) write(iout,9070) row,line,
     1                               (a(row,j),j=colmin,collim)
            endif
   90    continue
  100 continue
c
c
      return
      end
