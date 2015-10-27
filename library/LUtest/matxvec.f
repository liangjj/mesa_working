      subroutine matxvec(amat,colptr,rowind,nnz,fvec,svec,n2d)
c
c  multiplies the sparse  matrix amat times the vector fvec
c  result in vector svec.
c  matrix is indexed as for the linear equation solver SuperLU
c
c  n2d  = the order of the matrix 
c
      implicit complex*16(a-h,o-z)
      complex*16 fvec(n2d),svec(n2d)     
      complex*16 amat(nnz)
      integer colbeg, colend
      integer colptr(n2d+1),rowind(nnz),nnz
c
      do i=1,n2d
       svec(i) = (0.d0,0.d0)
      enddo
c
      do jcol=1,n2d
       colbeg = colptr(jcol)
       colend=colptr(jcol+1) -1
       do loop=colbeg,colend
           irow = rowind(loop)
           svec(irow)=svec(irow)+amat(loop)*fvec(jcol)
c
c to multiply by the TRANSPOSE of amat, replace the above statement by:
c
c           svec(jcol)=svec(jcol)+amat(loop)*fvec(irow)
c
       enddo
      enddo
c
      return
      end
