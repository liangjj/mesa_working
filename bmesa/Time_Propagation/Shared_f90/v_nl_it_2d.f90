!deck v_nl_it_2d.f
!**begin prologue     v_nl_it_2d
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           copy
!**author             schneider, barry (nsf)
!**source
!**purpose            mean time dependent hartree field
!**references
!**routines called
!**end prologue       v_nl_it_2d
  SUBROUTINE v_nl_it_2d
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE   
  REAL*8, DIMENSION(nphy(2),nphy(1))             :: wave_function
  REAL*8                                         :: fi, fj
  INTEGER                                        :: i, j, count
  count=0
  DO j=1,nphy(2)
     fj = grid(2)%f(j,j) * grid(2)%f(j,j)
     DO i=1,nphy(1)
        fi = grid(1)%f(i,i) * grid(1)%f(i,i)
        count=count+1
        v_tot(count) = v_tot(count)    +                                       &
                       nl_coef * wave_function(j,i) * wave_function(j,i)       &
                                                  *                            &
                                               fj * fi
     END DO
  END DO
  IF(log_main(5)) THEN
     title='non-linear potential'
     CALL prntrm(title,v_tot,n3d,1,n3d,1,iout)
  END IF
1 FORMAT(/,'calculating non-linear perturbing potential')
END SUBROUTINE v_nl_it_2d


