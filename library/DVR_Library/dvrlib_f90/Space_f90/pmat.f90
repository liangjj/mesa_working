! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Generalized First Derivative Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck pmat.f
!***begin prologue     pmat
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate generalized first derivative matrix elements.
!***references
!***routines called
!***end prologue       pmat
!\begin{eqnarray}
!\end{eqnarray}
  SUBROUTINE pmat(p_momr,far,dfbr,ptr,wtr,coord,nr,region)
  USE dvr_global,     ONLY   : parity, mass, iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nr
  REAL*8, DIMENSION(nr,nr)               :: p_momr, far, dfbr
  REAL*8, DIMENSION(nr)                  :: ptr, wtr
  CHARACTER (LEN=*)                      :: coord
  INTEGER                                :: region
  REAL*8                                 :: scale
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=3)                      :: itoc
  INTEGER                                :: i, j
!
  p_momr=0.d0
  IF(coord /= 'rho') THEN
     DO  i=1,nr
         p_momr(i,:) = p_momr(i,:) + far(i,i)*dfbr(i,:)*wtr(i)
     END DO
!
!    Add Bloch contributions
!
     p_momr(nr,nr) = p_momr(nr,nr) - far(nr,nr)*far(nr,nr)
     p_momr(1,1)   = p_momr(1,1)   + far(1,1)*far(1,1)
  ELSE
!     CALL lnkerr('not yet implimented')
  
!        In the cylindrical case, the functions are defined wrt the
!        argument $\rho$ or $\rho^{2}$ depending on the $m$ value
!        and we need to use the chain rule on the derivatives to get them
!        wrt $\rho$.  This is important in order to
!        get the limiting property at $\rho = 0$ correctly.
  
     IF(parity == 'even') THEN
!       need to put in code
     ELSE IF(parity == 'odd'.OR.parity == 'none') THEN
!       need to put in code      
     END IF
  END IF
  scale= .5D0
  p_momr = scale*p_momr
  IF(prn(3)) THEN
     title='unnormalized first derivative matrix for '// 'region = '//itoc(region)
     CALL prntrm(title,p_momr,nr,nr,nr,nr,iout)
  END IF
END SUBROUTINE pmat
