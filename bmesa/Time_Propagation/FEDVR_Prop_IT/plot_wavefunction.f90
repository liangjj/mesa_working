!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                       MODULE plot_wavefunction
!**begin prologue     plot_wavefunction
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       plot and print wavefunction
!
                         INTERFACE plot_psi
         MODULE PROCEDURE plot_psi_1d_d, plot_psi_2d_d, plot_psi_3d_d
                     END INTERFACE plot_psi             
                         INTERFACE print_psi
         MODULE PROCEDURE print_psi_1d_d, print_psi_2d_d, print_psi_3d_d
                     END INTERFACE print_psi             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         CONTAINS
!deck plot_psi_1d_d.f
!***begin prologue     plot_psi_1d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_psi_1d_d
  SUBROUTINE plot_psi_1d_d(vector,t)
  USE dvrprop_global_it
  USE dvr_shared  
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))             :: vector
  INTEGER                                :: t, i
  CHARACTER (LEN=16)                     :: fptoc
  IF (log_main(8)) THEN
      title='solution at t = '//fptoc(tim_pts(t+1))
      CALL print_psi(vector)
  END IF
  IF (plot) THEN
      DO i = 1,nphy(1)
         write(iplot(3),*) grid(1)%pt(i), vector(i)
      END DO
  END IF
END SUBROUTINE plot_psi_1d_d
!deck plot_psi_2d_d.f
!***begin prologue     plot_psi_2d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_psi_2d_d
  SUBROUTINE plot_psi_2d_d(vector,t)
  USE dvrprop_global_it
  USE dvr_shared  
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1))        :: vector
  INTEGER                                   :: t, i, j
  CHARACTER (LEN=16)                        :: fptoc
  IF (log_main(8)) THEN
      title='solution at t = '//fptoc(tim_pts(t+1))
      CALL print_psi(vector)
  END IF
  IF (plot) THEN
      DO i = 1,nphy(2)
         DO j = 1,nphy(1)
            write(iplot(3),*) grid(2)%pt(i), grid(1)%pt(j), vector(i,j)
         END DO
      END DO
  END IF
END SUBROUTINE plot_psi_2d_d
!deck plot_psi_3d_d.f
!***begin prologue     plot_psi_3d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_psi_3d_d
  SUBROUTINE plot_psi_3d_d(vector,t)
  USE dvrprop_global_it
  USE dvr_shared  
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))        :: vector
  INTEGER                                           :: i, j, k, t
  CHARACTER (LEN=16)                                :: fptoc
  IF (log_main(8)) THEN
      title='solution at t = '//fptoc(tim_pts(t+1))
      CALL print_psi(vector)
  END IF
  IF (plot) THEN
      DO i = 1,nphy(3)
         DO j = 1,nphy(2)
            DO k = 1,nphy(1)
               write(iplot(3),*) grid(3)%pt(i), grid(2)%pt(j),         &
                                 grid(3)%pt(k), vector(i,j,k)      
            END DO
         END DO
      END DO
   END IF
END SUBROUTINE plot_psi_3d_d
!deck print_psi_1d_d
!***begin prologue     print_psi_1d_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***references
!***routines called
!***end prologue       print_psi_1d_d
!
  SUBROUTINE print_psi_1d_d(wave_function)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: i
  REAL*8, DIMENSION(nphy(1))             :: wave_function
!
  write(iout,1) title
  write(iout,2)
  DO i=1,nphy(1)
     write(iout,3) grid(1)%pt(i), wave_function(i)
  END DO
1 FORMAT(a80)
2 FORMAT('     x      ',5x,'   Psi  ')
3 FORMAT(e12.6,5x,e12.6,5x)
END SUBROUTINE print_psi_1d_d
!deck print_psi_2d_d
!***begin prologue     print_psi_2d_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***references
!***routines called
!***end prologue       print_psi_2d_d
!
  SUBROUTINE print_psi_2d_d(wave_function)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                   :: i, j
  REAL*8, DIMENSION(nphy(2),nphy(1))        :: wave_function
!
  write(iout,1) title
  write(iout,2)
  DO i = 1,nphy(1)
     DO j = 1, nphy(2)
        write(iout,3) grid(1)%pt(i), grid(2)%pt(j),                 &
                      wave_function(j,i)
     END DO
  END DO
1 FORMAT(a80)
2 FORMAT('     x      ',5x,'     y      ',5x,'  Psi  ')
3 FORMAT(e12.6,5x,e12.6,5x,e12.6,5x)
END SUBROUTINE print_psi_2d_d
!deck print_psi_3d_d
!***begin prologue     print_psi_3d_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***references
!***routines called
!***end prologue       print_psi_3d_d
!
  SUBROUTINE print_psi_3d_d(wave_function)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                           :: i, j, k
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))        :: wave_function
!
  write(iout,1) title
  write(iout,2)
  DO i = 1,nphy(1)
     DO j = 1, nphy(2)
        DO k = 1, nphy(3)
           write(iout,3) grid(1)%pt(i), grid(2)%pt(j), grid(3)%pt(k),  &
                         wave_function(k,j,i)
        END DO
     END DO
  END DO
1 FORMAT(a80)
2 FORMAT('     x      ',5x,'     y      ',5x,'     z     ','         &
                                         Psi  ')
3 FORMAT(e12.6,5x,e12.6,5x,e12.6,5x,e12.6,5x,e12.6)
END SUBROUTINE print_psi_3d_d
END MODULE plot_wavefunction

