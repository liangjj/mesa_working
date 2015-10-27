! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Modify Diagonal DVR Matrix Elements for Time Propagation}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck modify_diag.f
!***begin prologue     modify_diag
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            modify diagonal DVR matrix elements
!***
!***references
!***routines called
!***end prologue       modify_diag
!
  SUBROUTINE modify_diag(ke,v,nf,modify,mat_typ)
  USE dvr_global
  USE dvrprop_global
  USE fd_global
  IMPLICIT NONE
  INTEGER                                :: nf, i
  REAL*8, DIMENSION(nf)                  :: v
  REAL*8, DIMENSION(nf,nf)               :: ke
  LOGICAL                                :: modify
  CHARACTER(LEN=*)                       :: mat_typ
!
!
  IF( mat_typ == 'full') THEN
     IF(modify) then
        DO i=1,nf
           ke(i,i) = ke(i,i) + v(i)
           v(i) = ke(i,i)
        END DO
        if(log_main(9)) then
            title='Modified Diagonal of Kinetic Energy Matrix. Potential Set To Zero'
            call prntfmn(title,v,nf,1,nf,1,iout,'e')
        END IF
        v = 0.d0
     ELSE
        DO i=1,nf
           v(i) = ke(i,i) + v(i)
        END DO
        if(log_main(9)) then
            title='Modified Potential Energy Matrix. Diagonal Kinetic Energy Set To Zero'
            call prntfmn(title,v,nf,1,nf,1,iout,'e')
        END IF
        DO i=1,nf
           ke(i,i) = 0.d0
        END DO
     END IF
  ELSE IF( mat_typ == 'banded') THEN
     IF(modify) THEN
        DO i=1,nf
           ke(i,1) = dscale*d(1) + v(i)
        END DO
        v = ke(:,1)
        if(log_main(9)) then
            title='Modified Diagonal of Kinetic Energy Matrix. Potential Set To Zero'
            call prntfmn(title,v,nf,1,nf,1,iout,'e')
        END IF
        v = 0.d0
     ELSE
        DO i=1,nf
           v(i) = v(i) + dscale*d(1)
        END DO
        if(log_main(9)) then
            title='Modified Potential Energy Matrix. Diagonal Kinetic Energy Set To Zero'
            call prntfmn(title,v,nf,1,nf,1,iout,'e')
        END IF
        ke(:,1) = 0.d0
     END IF
  ELSE
     call lnkerr('quit')
  END IF
  write(iplot(1),*) 'diagonal kinetic plus potential matrix'
  write(iplot(1),*) v
END SUBROUTINE modify_diag
