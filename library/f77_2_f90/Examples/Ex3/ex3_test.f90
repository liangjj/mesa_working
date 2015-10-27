  SUBROUTINE ex3_test(g,h,n,m)
  USE test_global,    ONLY    : offset
  real*8, dimension(n,m) :: g 
  real*8, dimension(n)   :: h 
  do i=1,n
   do j=1,m
     g(i,j)=i+j+offset
   END DO
   h(i)=i+offset            
  END DO
  RETURN
  END SUBROUTINE ex3_test
















