!deck v_nl_it_3d.f
!**begin prologue     v_nl_it_3d
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           copy
!**author             schneider, barry (nsf)
!**source
!**purpose            mean time dependent hartree field
!**references
!**routines called
!**end prologue       v_nl_it_3d
  SUBROUTINE v_nl_it_3d
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE   
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))             :: wave_function
  REAL*8                                                 :: fi, fj, fk
  INTEGER                                                :: i, j, k, count
  count=0
  DO k=1,nphy(3)
     fk = grid(3)%f(k,k) * grid(3)%f(k,k)
     DO j=1,nphy(2)
        fj = grid(2)%f(j,j) * grid(2)%f(j,j)
        DO i=1,nphy(1)
           fi = grid(1)%f(i,i) * grid(1)%f(i,i)
           count=count+1
           v_tot(count) = v_tot(count)    +                                           &
                          nl_coef * wave_function(k,j,i) * wave_function(k,j,i)       &
                                                       *                              &
                                                    fk * fj * fi
        END DO
     END DO
  END DO
  IF(log_main(5)) THEN
     title='non-linear potential'
     CALL prntrm(title,v_tot,n3d,1,n3d,1,iout)
  END IF
1 FORMAT(/,'calculating non-linear perturbing potential')
END SUBROUTINE v_nl_it_3d


