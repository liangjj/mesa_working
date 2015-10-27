! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Radial Kinetic Energy Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck kinrad.f
!***begin prologue     kinrad
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           kinetic energy
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements
!***                   in spherical coordinates.
!***
!***description
!***references
!***routines called
!***end prologue       kinrad
  SUBROUTINE kinrad(kmat,p,dp,ddp,pt,wt,n,parity)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n    
  REAL*8, DIMENSION(n,n)                 :: kmat, p, dp, ddp
  REAL*8, DIMENSION(n)                   :: pt, wt
  CHARACTER (LEN=*)                      :: parity
  REAL*8                                 :: two=2.d0, four=4.d0, six=6.d0
  INTEGER                                :: i, j 
  IF(parity == 'none') THEN
     DO  i=1,n
         kmat(i,:) = kmat(i,:) + p(i,i) * wt(i) *  &
                   ( ddp(i,:) + ( two/pt(i) ) * dp(i,:) )
     END DO
  ELSE IF(parity == 'even'.OR.parity == 'odd') THEN
     DO  i=1,n
         kmat(i,:) = kmat(i,:) + p(i,i) * wt(i) *  &
                   ( four * pt(i) * pt(i) * ddp(i,:) &
                         + six * dp(i,:) )
     END DO
  ELSE
     write(output,1)
     stop
  END IF
1 FORMAT('error in argument type')
END SUBROUTINE kinrad
