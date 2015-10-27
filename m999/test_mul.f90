subroutine test_mul(a,b,c,n)
      USE input_output
      implicit none
      INTEGER                                :: n, i, j
      REAL*8, DIMENSION(n,n)                 :: a, b, c
      REAL*4                                 :: secnds
      REAL*4                                 :: del_t
      REAL*4, DIMENSION(2)                   :: time
      time(1)=secnds(0.0)
      DO i=1,n
         DO j=1,i
            a(i,j) = 1.d0/(i+j-1)
            b(i,j) = a(i,j)
            a(j,i) = a(i,j)
            b(j,i) = a(j,i)
         END DO
      END DO
      time(2)=secnds(0.0)
      del_t=time(2)-time(1)
      write(iout,*) 'Time to Set Up Matrix = ', del_t
      time(1)=secnds(0.0)
      call ebc(c,a,b,n,n,n)
      time(2)=secnds(0.0)
      del_t=time(2)-time(1)
      write(iout,*) 'Time to Perform Matrix Mutiply = ', del_t
end subroutine

