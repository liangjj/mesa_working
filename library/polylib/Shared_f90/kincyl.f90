! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Cylindrical Kinetic Energy Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck kincyl.f
!***begin prologue     kincyl
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           kinetic energy
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements
!***                   in cartesian coordinates.
!***
!***description
!***references
!***routines called
!***end prologue       kincyl
  SUBROUTINE kincyl(kmat,p,dp,ddp,pt,wt,n,parity)
  USE dvr_global,   ONLY    : inp, iout
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n,n)                 :: kmat, p, dp, ddp
  REAL*8, DIMENSION(n)                   :: pt, wt
  CHARACTER (LEN=*)                      :: parity
  REAL*8                                 :: one=1.d0, two=2.d0, &
                                            four=4.d0, six=6.d0
  INTEGER                                :: i, j
  IF(parity == 'none') THEN
     DO  i=1,n
         kmat(i,:) = kmat(i,:)  &
                   + p(i,i) * wt(i) * ( ddp(i,:) +  &
                    (one/pt(i)) * dp(i,:) )
     END DO
  ELSE IF(parity == 'even'.OR.parity == 'odd') THEN
     DO  i=1,n
         kmat(i,:) = kmat(i,:)  &
                   + p(i,i) * wt(i) * four *  &
                   ( pt(i) * pt(i) * ddp(i,:) + dp(i,:) )
     END DO
  ELSE
     CALL lnkerr('error in argument type')
  END IF
END SUBROUTINE kincyl
