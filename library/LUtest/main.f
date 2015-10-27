c
c  Example for using SuperLU to solve linear equations -- McCurdy 2/19/00
c
c  Uses the standard sparse matrix indexing explained in
c  SuperLU documentation
c
c  The routiine "matrix" builds a matrix while computing
c  only the upper triangle, but storing, IN PACKED FORM, the
c  entire matrix.
c
c  The routine "matxvec" multiplies a packed matrix times a vector
c  and is used here to check the result.
c  
c  The routine "slu_fortran" is a driving routine that does three
c  things:  reorder the matrix for optimal bandwidth and factorization
c  properties,  factor and save factorization, solve.
c

      parameter(nmax=200,nonzmax=nmax*100)
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)
c
c dimension arrays for SuperLU
c
      complex*16 zlineq(nonzmax)
      integer*4 colptr(nmax+1),rowind(nonzmax)
      integer*4 perm_c(nmax)
      complex*16  psi(nmax),psisave(nmax),psicomp(nmax)
c
c  n2d is the order of the matrix n2d < nmax
c
      n2d =15 
c
c  construct (and save)  right hand side
c
      do i=1, n2d
       psi(i) = i
       psisave(i) = psi(i)
      enddo
c
c  construct the matrix 
c
      call matrix(zlineq,n2d,nonzmax,colptr,rowind,nnz)
      write(6,101) n2d, nnz
  101 format(//,'Matrix of order',i5,' with ',i8,' nonzeros',//)
c
c
c  solve linear equations for vector psi  
c
c  call to SuperLU driving routine "SLU_fortran" written by
c  Mark Baertschy
c
      nrhs = 1
      ldb = n2d
       call SLU_fortran(n2d,nnz,nrhs,zlineq,rowind,colptr,perm_c, 
     $  psi,ldb,info,'Y')
c
c if there was in fault in SuperLU, write a warning 
c
      if(info.ne.0) then
        write(6,122) info
  122   format('crash in SuperLU (linear solver) ',i10)
       endif
c
c next call is just to free up memory
c
       call SLU_fortran(n2d,nnz,nrhs,zlineq,rowind,colptr,perm_c, 
     $  psi,ldb,info,'F')
c
      write(6,200) 
  200 format(//,'   Solution vector ')
      write(6,201) (psi(iii),iii=1,n2d)
  201 format(1x,2e12.5)
c 
c  check solution
c
      call matrix(zlineq,n2d,nonzmax,colptr,rowind,nnz)
      write(6,202) 
  202 format(//,' matrix times  solution                rhs') 
      call matxvec(zlineq,colptr,rowind,nnz,psi,psicomp,n2d)
      do i=1,n2d
       write(6,104) psicomp(i), psisave(i)
 104  format(2e12.5,5x,2e12.5)
      enddo
c
      stop 
      end    
