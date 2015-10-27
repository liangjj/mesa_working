! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{SOLN: Soln Vectors}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck soln
!**begin prologue     soln
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            combine initial solution with solution to driven
!**                   equation to get total final solution and write to disk.
!**references
!***routines called
!***end prologue       soln
  SUBROUTINE soln(tim,rtemp)
  USE arnoldi_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER                               :: tim, i
  CHARACTER (LEN=16)                    :: fptoc
  REAL*8, DIMENSION(n3d)                :: rtemp 
  COMPLEX*16                            :: conjg       
  chi = chi + psi0
  IF(log_main(8)) then
     title='solution at t ='//fptoc(t1)
     CALL prntcmn(title,chi,n3d,1,n3d,1,iout,'e')
  END IF
  CALL iosys ('write real solution to bec',n3d*2,chi,0,' ')
!  CALL iosys('write real wavefunction to plot '//  &
!             'without rewinding',2*n3d,chi,0,' ')
  IF(plot.and.tim == ntreg) then
     do i=1,n3d
        write(iplot(3),*) grid(1)%pt(i), real(chi(i))
        write(iplot(4),*) grid(1)%pt(i), imag(chi(i))
        rtemp(i) =  sqrt( chi(i)*conjg(chi(i)) )
        write(iplot(5),*) grid(1)%pt(i), rtemp(i)
     END DO
  END IF
  call auto(auto_corr(tim),soln_0,chi)
1 FORMAT(a16)
2 FORMAT(6e15.8)
END SUBROUTINE soln





