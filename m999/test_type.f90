program test_type
      USE input_output
      implicit none
      REAL*8, dimension(:,:), allocatable    :: a, b, c
      TYPE mat
         REAL*8, DIMENSION(:,:), POINTER     :: ma, mb, mc
      END TYPE mat
      TYPE(mat)                              :: mat_pnt
      INTEGER                                :: intkey, n_size
      CHARACTER(LEN=4096)                    :: ops
      call drum
      call IOsys('read character options from rwf',-1,0,0,ops)
      n_size=intkey(ops,'matrix-size',100,' ') 
      write(iout,*) 'Matrix Size = ',n_size

!      ALLOCATE(a(n_size,n_size),b(n_size,n_size),c(n_size,n_size))
      ALLOCATE(mat_pnt%ma(n_size,n_size),mat_pnt%mb(n_size,n_size),mat_pnt%mc(n_size,n_size))
      CALL test_mul(mat_pnt%ma,mat_pnt%mb,mat_pnt%mc,n_size)
!      CALL test_mul(a,b,c,n_size)
      DEALLOCATE(mat_pnt%ma,mat_pnt%mb,mat_pnt%mc)
!      DEALLOCATE(a,b,c)
      call chainx(0)
      stop  
end program test_type
