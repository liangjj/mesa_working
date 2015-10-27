      subroutine matrix(amat,n2d,nonzmax,colptr,rowind,nnz)
c
c  n is the order of the matrix
c  band sets the bandwidth
c
c  nonzmax is maximum number of nonzero elements
c
      implicit complex*16 (a-h,o-z)
      integer band
      dimension amat(nonzmax)
      integer colptr(n2d+1), colbeg,colend
      integer rowind(nonzmax)
c
      band = 4
c
c
c initialize count of nonzero matrix elements
      icount = 0
c
c loop on first index of matrix 
c
      do 50 i2d=1,n2d
      colptr(i2d) = icount+1
c
c loop on second index of matrix.  matrix  is computed
c in this loop
c
      do 49 j2d=1,n2d
c
c
c check if we are in the upper triangle, if so,
c find previously computed matrix element 
c
      if(j2d.lt.i2d) then
       colbeg = colptr(j2d)
       colend = colptr(j2d+1)-1
         do irow=colbeg,colend
           index=rowind(irow)
           if(index.eq.i2d) then
              val = amat(irow)
              icount=icount+1
              amat(icount) = val
              rowind(icount) = j2d
           endif
          enddo
c skip out of j2d loop after finding the right matrix element
c previously computed if we are in the upper triangle
         go to 49
        endif
c
c  ask if the element is nonzero
c
          if (iabs(i2d-j2d).gt.band) go to 49
c
c   if not, "compute" its value
c
          val =  dcmplx(float(i2d),float(j2d)) 
          icount=icount+1
c
           if(icount.gt.nonzmax) then
           write(6,999) icount, nonzmax
  999      format('OUCH! more elements than nonzmax!',2i14)
           stop
           endif
c
c  load the value in the array, using sparse indexing
c  and load the corresponding row index in rowind
c          
           amat(icount) = val
           rowind(icount) = j2d
c
  49  continue
  50  continue
c
c finish the column pointer array and set the count of
c nonzero elements
c
           colptr(n2d+1) = icount+1
           nnz = icount
c
       return
       end

