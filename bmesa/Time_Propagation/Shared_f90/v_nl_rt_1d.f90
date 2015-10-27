!deck v_nl_rt_1d.f
!**begin prologue     v_nl_rt_1d
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           copy
!**author             schneider, barry (nsf)
!**source
!**purpose            mean time dependent hartree field
!**references
!**routines called
!**end prologue       v_nl_rt_1d
  SUBROUTINE v_nl_rt_1d(wave_function)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE   
  REAL*8, DIMENSION(nphy(1),2)           :: wave_function
  INTEGER                                :: i
  DO i=1,nphy(1)
     v_tot(i) = v_tot(i) + nl_coef * (wave_function(i,1) * wave_function(i,1)       &
                                                       +                            &
                                    wave_function(i,2) * wave_function(i,2) )       &
                                                       *                            &
                                        grid(1)%f(i,i) * grid(1)%f(i,i)
  END DO
  IF(log_main(5)) THEN
     title='non-linear potential'
     CALL prntrm(title,v_tot,n3d,1,n3d,1,iout)
  END IF
1 FORMAT(/,'calculating non-linear perturbing potential')
END SUBROUTINE v_nl_rt_1d


