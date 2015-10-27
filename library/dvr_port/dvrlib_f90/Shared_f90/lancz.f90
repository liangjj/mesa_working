! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Lanczos}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck lancz.f
!***begin prologue     lancz
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           eigenvalues, eigenvectors
!***author             schneider, barry (nsf)
!***source
!***purpose            generate orthogonal polynomials using recursion.
!***                   using the generated alpha and beta coefficients find
!***                   the points and weights of the generalized gauss
!***                   quadrature.
!***
!***description
!***references
!***routines called
!***end prologue       lancz
!\begin{eqnarray}
!  \beta_{j} P_{j}(x) &=& \big ( x - \alpha_{j} \big ) P_{j-1}(x) - \beta_{j-1} P_{j-2}(x)
!                             \\ \nonumber
!  \alpha_{j} &=& \langle P_{j-1} \mid x \mid P_{j-1} \rangle \\ \nonumber
!  \beta_{j} &=& \langle P_{j-1} \mid x \mid P_{j} \rangle \\ \nonumber
!  \beta_{j} P_{j}^{\prime}(x) &=& \big ( x - \alpha_{j} \big ) P_{j-1}^{\prime}(x)
!                             - \beta_{j-1} P_{j-2}^{\prime}(x) + P_{j-1}(x)
!                                 \\ \nonumber
!  \beta_{j} P_{j}^{\prime \prime}(x) &=& \big ( x - \alpha_{j} \big )
!                                      P_{j-1}^{\prime \prime}(x)
!                            - \beta_{j-1} P_{j-2}^{\prime \prime}(x)
!                            + 2 P_{j-1}^{\prime}(x)
!\end{eqnarray}
  SUBROUTINE lancz(v,arg,a,b,wt,refwt,scr,n,iter)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n, iter
  REAL*8, DIMENSION(n,0:iter)            :: v
  REAL*8, DIMENSION(n)                   :: arg, wt, refwt, scr
  REAL*8, DIMENSION(iter)                :: a, b
  REAL*8                                 :: anorm
  INTEGER                                :: i
  REAL*8, DIMENSION(:), ALLOCATABLE      :: vtmp
  ALLOCATE(vtmp(n))
!     first vector is input.  normalize.
  vtmp=v(:,0)*v(:,0)*wt*refwt
  anorm=SQRT(1.d0/sum(vtmp,1))
  v(:,0)=anorm*v(:,0)
!     we now have the first function
  IF (iter > 1) THEN
!         form argument times the first function
      scr=arg*v(:,0)
!         calculate a(1)
      vtmp=v(:,0)*scr*wt*refwt
      a(1)=sum(vtmp,1)
!         form mat times the first function - a(1) times the first function
!         and store it in the next polynomial
      v(:,1)=scr - a(1)*v(:,0)
!         calculate b(1)
      vtmp=v(:,1)*v(:,1)*wt*refwt
      b(1)=SQRT( sum(vtmp,1) )
!         normalize the second polynomial
      v(:,1)=(1.d0/b(1))*v(:,1)
  END IF
  IF (iter > 2) THEN
      DO  i=2,iter
!            multiply the last calculated polynomial by mat
          scr=arg*v(:,i-1)
!            orthogonalize it to the two previous polynomials
!            calculating a(i) as we go
          vtmp=v(:,i-1)*scr*wt*refwt
          a(i)=sum(vtmp,1)
          v(:,i)=scr - a(i)*v(:,i-1) - b(i-1)*v(:,i-2)
!            calculate b(i)
          vtmp=v(:,i)*v(:,i)*wt*refwt       
          b(i) = SQRT(sum(vtmp,1))
!            normalize the polynomial and we are done
          v(:,i)=(1.d0/b(i))*v(:,i)
      END DO
  END IF
  DEALLOCATE(vtmp)
END SUBROUTINE lancz


