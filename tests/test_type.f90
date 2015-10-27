program test_type
      implicit none
      REAL*8, dimension(:,:), allocatable    :: a, b, c
      TYPE mat
         REAL*8, DIMENSION(:,:), POINTER     :: ma, mb, mc
      END TYPE mat
      TYPE(mat)                              :: mat_pnt
      INTEGER                                :: n_size
      n_size= 1000000
      ALLOCATE(a(n_size,n_size),b(n_size,n_size),c(n_size,n_size))
      CALL test_mul(a,b,c,n_size)
      DEALLOCATE(a,b,c)
end program test_type
subroutine test_mul(a,b,c,n)
      implicit none
      INTEGER                                :: n
      REAL*8, DIMENSION(n,n)                 :: a, b, c
end subroutine

