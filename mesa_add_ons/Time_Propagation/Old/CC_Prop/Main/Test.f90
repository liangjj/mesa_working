!
!deck Test
!**begin prologue     
!**date written       070420   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source
!**purpose            test
!**description        
!**                   
!**                   
!**references
!**routines called        
!                     
!
!**end prologue       Test
  PROGRAM Test
!
  IMPLICIT NONE
  INTEGER                                  :: inp, out, size, i, j, k
  REAL*8, DIMENSION(:,:), ALLOCATABLE      :: a, b, c
!
! Read in the Input
  inp=5
  out=6
  OPEN(inp,file='inp_test',status='old')
  OPEN(out,file='out_test',status='unknown')
!
  read(inp,*) size
  ALLOCATE ( a(size,size), b(size,size), c(size,size) )
  DO i=1,size
     DO j=1,i
        a(i,j) = i + j
        a(j,i) = a(i,j)
     END DO
  END DO
  b(:,:) = a(:,:)
  c(:,:) = 0.d0
  DO i=1,size
     DO j=1,size
        DO k=1,size
           c(i,j) = c(i,j) + a(i,k) * b(k,j)
        END DO
     END DO
  END DO
  call dgemm('n','n',size,size,size,1.d0,a,size,b,size,0.d0,c,size)
  write(out,1) c
1 Format(/,5e15.8)
  stop
END PROGRAM Test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
