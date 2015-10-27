! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Sc_dip}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck sc_dip.f
!**begin prologue     sc_dip
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           copy
!**author             schneider, barry (nsf)
!**source
!**purpose            compute spatial factors for dipole potentials.
!**references
!**routines called
!**end prologue       sc_dip
  SUBROUTINE sc_dip(tfac)
  USE Iterative_Global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8                                 :: tfac,  fac
  INTEGER                                :: i, j, k, count
  REAL*8, DIMENSION(:), ALLOCATABLE      :: scr
  IF(spdim == 1) THEN
     v_tot = v_tot + tfac*grid(1)%pt
  ELSE IF(spdim == 2) THEN
     ALLOCATE(scr(nphy(1)*nphy(2)))
     count=0
     DO  i=1,nphy(1)
         DO  j=1,nphy(2)
             count =count + 1
             scr(count) = SQRT(grid(1)%pt(i)*grid(1)%pt(i) + &
                           grid(2)%pt(j)*grid(2)%pt(j))
         END DO
     END DO
     v_tot = v_tot + tfac*scr   
     DEALLOCATE(scr)
  ELSE IF(spdim == 3) THEN
     ALLOCATE(scr(nphy(1)*nphy(2)*nphy(3)))
     count=0
     DO  i=1,nphy(1)
         DO  j=1,nphy(2)
             DO  k=1,nphy(3)
                 count = count + 1 
                 scr(count) = SQRT(grid(1)%pt(i)*grid(1)%pt(i) + &
                              grid(2)%pt(j)*grid(2)%pt(j) + &
                              grid(3)%pt(k)*grid(3)%pt(k))
             END DO
         END DO
     END DO
     v_tot = v_tot + tfac*scr   
     DEALLOCATE(scr)
  ELSE
    write(iout,1)
    stop
  END IF
1 format(/,5x,'error in sc_dip')
END SUBROUTINE sc_dip


