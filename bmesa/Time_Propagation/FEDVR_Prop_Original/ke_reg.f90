! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Make Diagonal DVR Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck ke_reg.f
!***begin prologue     ke_reg
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            regional matrix elements of kinetic energy           
!***                   and potential matrix elements
!***references
!***routines called
!***end prologue       ke_reg
!
  SUBROUTINE ke_reg(k_mat_d,n_tot,idim,mat_typ)
  USE dvr_global
  USE dvr_shared
  USE dvrprop_global
  IMPLICIT NONE
  INTEGER                                :: n_tot, idim
  REAL*8, DIMENSION(n_tot,n_tot)         :: k_mat_d
  INTEGER                                :: start, reg
  CHARACTER(LEN=*)                       :: mat_typ
!
!
! Calculate the needed matrix elements
  IF(mat_typ=='full') THEN
     start=1
     DO reg=1,nreg
        CALL fil_reg_dvr(k_mat_d(start,start),         &
                         mat_reg_d(reg,idim)%ke_mat_d, &
                         nfun_reg(reg,idim),nphy(idim),reg)
        start = start + nfun_reg(reg,idim) - 1
     END DO
  ELSE IF(mat_typ=='banded') THEN
     start=1
     DO reg=1,nreg
        call fil_reg_fd(k_mat_d(start,1),              &
                        mat_reg_d(reg,idim)%ke_mat_d,  &
                        nfun_reg(reg,idim),reg)
        start = start + nfun_reg(reg,idim) - 1
     END DO
  ELSE
     CALL lnkerr('error')
  END IF  
END SUBROUTINE ke_reg
