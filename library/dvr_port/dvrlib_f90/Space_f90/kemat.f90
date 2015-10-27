! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Generalized Kinetic Energy Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck kemat.f
!***begin prologue     kemat
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements with
!***                   singularities removed and derivative bloch operator
!***                   added.
!***
!***description
!***references
!***routines called
!***end prologue       kemat
!\begin{eqnarray}
!\end{eqnarray}
  SUBROUTINE kemat(ker,far,dfbr,ddfbr,ptr,wtr,coord,nr,region)
  USE dvr_global,     ONLY   : parity, mass, output
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nr
  REAL*8, DIMENSION(nr,nr)               :: ker, far, dfbr, ddfbr
  REAL*8, DIMENSION(nr)                  :: ptr, wtr
  CHARACTER (LEN=*)                      :: coord
  INTEGER                                :: region
  REAL*8                                 :: scale
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=3)                      :: itoc
  INTEGER                                :: i, j
!
  ker=0.d0
  IF(coord /= 'rho') THEN
     DO  i=1,nr
         ker(i,:) = ker(i,:) + far(i,i)*ddfbr(i,:)*wtr(i)
     END DO
!
!    Add Bloch contributions
!
     ker(nr,:) = ker(nr,:) - far(nr,nr)*dfbr(nr,:)
     ker(1,:)  = ker(1,:)  + far(1,1)*dfbr(1,:)
  ELSE
  
!        In the cylindrical case, the functions are defined wrt the
!        argument $\rho$ or $\rho^{2}$ depending on the $m$ value
!        and we need to use the chain rule on the derivatives to get them
!        wrt $\rho$.  This is important in order to
!        get the limiting property at $\rho = 0$ correctly.
  
     IF(parity == 'even') THEN
        DO  i=1,nr
            ker(i,:) = ker(i,:)  + far(i,i)*wtr(i)*4.d0*  &
                          ( ptr(i)*ptr(i)*ddfbr(i,:) + dfbr(i,:) )
        END DO
        ker(nr,:) = ker(nr,:) - 2.d0*ptr(nr)*ptr(nr)*far(nr,nr)*dfbr(nr,:)
        ker(1,:)  = ker(1,:)  + 2.d0*ptr(1)*ptr(1)*far(1,1)*dfbr(1,:)      
     ELSE IF(parity == 'odd'.OR.parity == 'none') THEN
        DO  i=1,nr
            ker(i,:) = ker(i,:)  + far(i,i)*wtr(i)* ( ddfbr(i,:) &
                                 + dfbr(i,:)/ptr(i) )
        END DO
        ker(nr,:) = ker(nr,:) - ptr(nr)*far(nr,nr)*dfbr(nr,:)
        ker(1,:)  = ker(1,:)  + ptr(1)*far(1,1)*dfbr(1,:)
     END IF
  END IF
  scale= -.5D0/mass
  ker=scale*ker
  IF(prn(3)) THEN
     title='unnormalized kinetic energy matrix for '// 'region = '//itoc(region)
     CALL prntrm(title,ker,nr,nr,nr,nr,output)
  END IF
END SUBROUTINE kemat
