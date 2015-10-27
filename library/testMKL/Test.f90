Program Test
IMPLICIT NONE
REAL*8, DIMENSION(100)  :: a
REAL*8, DIMENSION(100)  :: b
REAL*8                  :: ddot
REAL*8                  :: prod
a(:) = 1.d0
b(:) = .5d0
prod=ddot(100,a,1,b,1)
write(6,*) 'Dotpro = ',prod
END PROGRAM Test
