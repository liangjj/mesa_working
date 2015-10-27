! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Kinetic Energy Matrix Elements with Hermite Weight Function}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck ke_hermite.f
!***begin prologue     ke_hermite
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           kinetic energy with Hermite weight function
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements
!***                   for a hermite weight function.
!***
!***description
!***references
!***routines called
!***end prologue       ke_hermite
  SUBROUTINE ke_hermite(ke,fa,fb,dfb,ddfb,pt,wt,n,region)
  USE dvr_global,     ONLY   : mass, output
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: n, region    
  REAL*8, DIMENSION(n,n)                 :: ke, fa, fb
  REAL*8, DIMENSION(n,n)                 :: dfb, ddfb
  REAL*8, DIMENSION(n)                   :: pt, wt
  REAL*8                                 :: scale
  REAL*8                                 :: one=1.d0, two=2.d0
  INTEGER                                :: i 
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=3)                      :: itoc
  DO  i=1,n
      ke(i,:) = ke(i,:) + fa(i,i) * wt(i) * &
                  ( ddfb(i,:) - two * pt(i) * dfb(i,:) )
  END DO
  DO i=1,n
     ke(i,i) = ke(i,i) + fa(i,i) * wt(i) * ( pt(i)*pt(i) - one ) * fb(i,i)
  END DO
  scale= -.5D0/mass
  ke=scale*ke
  IF(prn(3)) THEN
     title='unnormalized kinetic energy matrix for ' &
            // 'region = '//itoc(region)
     CALL prntrm(title,ke,n,n,n,n,output)
  END IF
END SUBROUTINE ke_hermite


