! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Construct Final DVR Potential}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck conv.f
!***begin prologue     conv
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            DVR/FEM basis.
!***
!***references

!***routines called
!***end prologue       conv

  SUBROUTINE conv(vmat,mat,v,nglobal)
!
  USE dvr_global,    ONLY   : nreg, npt, iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                            :: nglobal
  INTEGER, DIMENSION(nreg)           :: v               
  REAL*8, DIMENSION(nglobal)         :: vmat
  REAL*8, DIMENSION(*)               :: mat
  CHARACTER (LEN=80)                 :: title
  INTEGER                            :: reg, start, end, nfun, last, row 
!
  vmat=0.d0
  row=1
  DO  reg=1,nreg
      start=2
      END=npt(reg)
      IF(reg == 1) THEN
         start=1
      END IF
      nfun = END - start + 1
      last = row + nfun - 1
      CALL filv(vmat(row),mat(v(reg)),npt(reg),nfun, &
                nglobal,start)
      row = row + nfun
  END DO
  IF(prn(1)) THEN
     title='global potential'
     CALL prntrm(title,vmat,nglobal,1,nglobal,1,iout)
  END IF
1    FORMAT(/,5X,'region = ',i3,1X,'number of functions = ',i5,/,5X,  &
    'starting function = ',i5,1X,'ending function = ',i5,  &
    /,5X,'global starting function = ',i5,1X, 'global ending function = ',i5,  &
    /,5X,'bridge function = ',l1)
END SUBROUTINE conv




