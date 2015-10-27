program testblas
      implicit none
      real*8, dimension(:), allocatable    :: a_vec, b_vec
      real*8, dimension(:,:), allocatable  :: a_mat, b_mat, c_mat
      integer                              :: n, i, j
      real*8                               :: sdot
      real*8                               :: result
      n=50
      ALLOCATE(a_vec(n),b_vec(n))
      a_vec=1.d0
      b_vec=2.d0
      result=sdot(n,a_vec,1,b_vec,1)
      write(6,1) result
      DEALLOCATE(a_vec,b_vec)
      ALLOCATE(a_mat(n,n),b_mat(n,n),c_mat(n,n))     
      a_mat=0.d0
      b_mat=0.d0
      c_mat=0.d0
      DO i=1,n
         DO j=1,n
            b_mat(i,j)=1.0d0
            c_mat(i,j)=1.d0
         END DO
      END DO
      call ebc(a_mat,b_mat,c_mat,n,n,n)
      write(6,*) a_mat
      DEALLOCATE(a_mat,b_mat,c_mat)
1 Format('Result = ',e15.8)
end program testblas
