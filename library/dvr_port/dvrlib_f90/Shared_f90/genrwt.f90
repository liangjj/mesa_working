! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Reference Weight Function}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck genrwt.f
!***begin prologue     genrwt
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           reference weight
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate weight functions and derivatives.
!***
!***description
!***
!***
!***references
!***routines called
!***end prologue       genrwt
  SUBROUTINE genrwt(rwt,drwt,ddrwt,pt,wt,alpha,beta,deriv,edge,n)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: rwt, drwt, ddrwt, pt
  CHARACTER (LEN=*)                      :: wt
  REAL*8                                 :: alpha, beta
  LOGICAL                                :: deriv
  REAL*8, DIMENSION(2)                   :: edge
  REAL*8, DIMENSION(:), ALLOCATABLE      :: fac1, fac2
!     Weight functions and their derivatives are computed for a variety of
!     cases
  IF(wt == 'one'.or.wt == 'legendre') THEN
     rwt = 1.d0
     IF(deriv) THEN
        drwt=0.d0
        ddrwt=0.d0
     END IF
  ELSE IF(wt == 'r') THEN
     rwt=pt
     IF(deriv) THEN
        drwt=1.d0
        ddrwt=0.d0
     END IF
  ELSE IF(wt == 'rr') THEN
    rwt=pt*pt 
    IF(deriv) THEN
       drwt=2.d0*pt
       ddrwt=2.d0
    END IF
  ELSE IF(wt == 'hermite') THEN
    rwt=EXP(-pt*pt)
    IF(deriv) THEN
       drwt = -2.d0*rwt*pt
       ddrwt= ( -2.d0 + 4.d0*pt*pt ) * rwt
    END IF
  ELSE IF(wt == 'chebyshev-1') THEN
    rwt= SQRT ( 1.d0/( 1.d0 - pt*pt ) )
    IF(deriv) THEN
       ddrwt = rwt*rwt*rwt
       drwt = pt*ddrwt
       ddrwt = ddrwt +  3.d0*pt*pt*ddrwt*rwt*rwt
    END IF
  ELSE IF(wt == 'chebyshev-2') THEN
    rwt= SQRT (  1.d0 -pt*pt )
    IF(deriv) THEN
       drwt = - pt/rwt
       ddrwt = ( - 1.d0 + pt*drwt/rwt )/rwt
    END IF
  ELSE IF(wt == 'laguerre') THEN
    rwt= pt**alpha*EXP(-pt)
    IF(deriv) THEN
       ALLOCATE(fac1(n),fac2(n))
       fac1 =  - pt + alpha/pt 
       fac2 =  1.d0 + alpha/(pt*pt) 
       drwt = fac1 * rwt
       ddrwt = fac1 * drwt - fac2 * rwt
       DEALLOCATE(fac1,fac2)
    END IF
  ELSE IF(wt == 'jacobi') THEN
    rwt= (1.d0-pt)**alpha * (1.d0+pt)**beta
    IF(deriv) THEN
       ALLOCATE(fac1(n),fac2(n))
       fac1= -alpha/(1.d0-pt) 
       fac2=beta/(1.d0+pt) 
       drwt = (fac1 + fac2 )*rwt
       ddrwt = ( fac1 + fac2 )*drwt
       fac1=fac1/(1.d0-pt)
       fac2=fac2/(1.d0+pt)
       ddrwt = ddrwt + ( fac1 + fac2 )*rwt
       DEALLOCATE(fac1,fac2)
    END IF
  ELSE IF(wt == 'rys') THEN
    rwt=EXP(-alpha*pt*pt)
    IF(deriv) THEN
      ALLOCATE(fac1(n))
      fac1 = -2.d0*pt*alpha
      drwt = fac1*rwt
      ddrwt = fac1*drwt - 2.d0*alpha*rwt
      DEALLOCATE(fac1)
    END IF
  ELSE
    rwt=1.d0
  END IF
END SUBROUTINE genrwt










