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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        INTERFACE plot_psi
             MODULE PROCEDURE plot_psi_1d_z, plot_psi_2d_z, plot_psi_3d_z
                    END INTERFACE plot_psi
                        INTERFACE print_psi
             MODULE PROCEDURE print_psi_1d_z, print_psi_2d_z,    &
                              print_psi_3d_z
                    END INTERFACE print_psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck plot_psi_1d_z.f
!***begin prologue     plot_psi_1d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_psi_1d_z
  SUBROUTINE plot_psi_1d_z(vector,t)
  USE dvrprop_global_rt
  USE dvr_shared  
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: vector
  REAL*8                                 :: rtemp
  INTEGER                                :: t, i, iostat
  CHARACTER (LEN=16)                     :: fptoc
  IF (log_main(8)) THEN
      title='solution at t = '//fptoc(tim_pts(t+1))
      CALL print_psi(vector)
  END IF
!  IF (plot.and.t == ntreg) THEN
  IF (plot) THEN
      OPEN (UNIT=iplot(3),FILE='real_wavefunction-'//fptoc(tim_pts(t+1)),       &
            ACCESS='sequential',FORM='formatted',                               &
            IOSTAT=IOSTAT,STATUS='unknown')
     IF(IOSTAT /= 0) THEN
        CALL lnkerr('error in file handling')
     END IF
     OPEN (UNIT=iplot(4),FILE='imaginary_wavefunction-'//fptoc(tim_pts(t+1)),   &
           ACCESS='sequential',FORM='formatted',         &
           IOSTAT=IOSTAT,STATUS='unknown')
     IF(IOSTAT /= 0) THEN
        CALL lnkerr('error in file handling')
     END IF
     OPEN (UNIT=iplot(5),FILE='absolute_wavefunction-'//fptoc(tim_pts(t+1)),   &
           ACCESS='sequential',FORM='formatted',                               &
           IOSTAT=IOSTAT,STATUS='unknown')
     IF(IOSTAT /= 0) THEN
        CALL lnkerr('error in file handling')
     END IF
      DO i = 1,nphy(1)
         rtemp = SQRT ( ( vector(i,1) * vector(i,1)        &
                                      +                    &
                          vector(i,2) * vector(i,2) ) )
         write(iplot(3),*) grid(1)%pt(i), vector(i,1)
         write(iplot(4),*) grid(1)%pt(i), vector(i,2)
         write(iplot(5),*) grid(1)%pt(i), rtemp
      END DO

  END IF
 write(iplot(3),*) 
 write(iplot(4),*) 
 write(iplot(5),*) 
END SUBROUTINE plot_psi_1d_z
!deck plot_psi_2d_z.f
!***begin prologue     plot_psi_2d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_psi_2d_z
  SUBROUTINE plot_psi_2d_z(vector,t)
  USE dvrprop_global_rt
  USE dvr_shared  
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),2)   :: vector
  REAL*8                                 :: rtemp
  INTEGER                                :: t, i, j
  CHARACTER (LEN=16)                     :: fptoc
  IF (log_main(8)) THEN
      title='solution at t = '//fptoc(tim_pts(t+1))
      CALL print_psi(vector)
  END IF
  IF (plot.and.t == ntreg) THEN
      DO i = 1,nphy(2)
         DO j = 1,nphy(1)
            rtemp = SQRT ( ( vector(i,j,1) * vector(i,j,1)        &
                                           +                      &
                             vector(i,j,2) * vector(i,j,2) ) )
            write(iplot(3),*) grid(2)%pt(i), grid(1)%pt(j), vector(i,j,1)
            write(iplot(4),*) grid(2)%pt(i), grid(1)%pt(j), vector(i,j,2)
            write(iplot(5),*) grid(2)%pt(i), grid(1)%pt(j), rtemp
         END DO
      END DO
  END IF
END SUBROUTINE plot_psi_2d_z
!deck plot_psi_3d_z.f
!***begin prologue     plot_psi_3d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_psi_3d_z
  SUBROUTINE plot_psi_3d_z(vector,t)
  USE dvrprop_global_rt
  USE dvr_shared  
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: vector
  REAL*8                                         :: rtemp
  INTEGER                                        :: i, j, k, t
  CHARACTER (LEN=16)                             :: fptoc
  IF (log_main(8)) THEN
      title='solution at t = '//fptoc(tim_pts(t+1))
      CALL print_psi(vector)
  END IF
  IF (plot.and.t == ntreg) THEN
      DO i = 1,nphy(3)
         DO j = 1,nphy(2)
            DO k = 1,nphy(1)
               rtemp = SQRT ( ( vector(i,j,k,1) * vector(i,j,k,1)      &
                                                +                      &
                                vector(i,j,k,2) * vector(i,j,k,2) ) )
               write(iplot(3),*) grid(3)%pt(i), grid(2)%pt(j),         &
                                 grid(1)%pt(k), vector(i,j,k,1)      
               write(iplot(4),*) grid(3)%pt(i), grid(2)%pt(j),         &
                                 grid(3)%pt(k), vector(i,j,k,2)
               write(iplot(5),*) grid(3)%pt(i), grid(2)%pt(j),         &
                                 grid(3)%pt(k), rtemp
            END DO
         END DO
      END DO
   END IF
END SUBROUTINE plot_psi_3d_z
!deck print_psi_1d_z
!***begin prologue     print_psi_1d_z
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***references
!***routines called
!***end prologue       print_psi_1d_z
!
  SUBROUTINE print_psi_1d_z(wave_function)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: i
  REAL*8, DIMENSION(nphy(1),2)           :: wave_function
!
  write(iout,1) title
  write(iout,2)
  DO i=1,nphy(1)
     write(iout,3) grid(1)%pt(i), wave_function(i,1), wave_function(i,2)
  END DO
1 FORMAT(a80)
2 FORMAT('     x      ',5x,'   Real Psi  ',5x,'Imag Psi ')
3 FORMAT(e12.6,5x,e12.6,5x,e12.6)
END SUBROUTINE print_psi_1d_z
!deck print_psi_2d_z
!***begin prologue     print_psi_2d_z
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***references
!***routines called
!***end prologue       print_psi_2d_z
!
  SUBROUTINE print_psi_2d_z(wave_function)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: i, j
  REAL*8, DIMENSION(nphy(2),nphy(1),2)   :: wave_function
!
  write(iout,1) title
  write(iout,2)
  DO i = 1,nphy(1)
     DO j = 1, nphy(2)
        write(iout,3) grid(1)%pt(i), grid(2)%pt(j),                 &
                      wave_function(j,i,1), wave_function(j,i,2)
     END DO
  END DO
1 FORMAT(a80)
2 FORMAT('     x      ',5x,'     y      ',5x,'  Real Psi  ',5x,' Imag Psi ')
3 FORMAT(e12.6,5x,e12.6,5x,e12.6,5x,e12.6)
END SUBROUTINE print_psi_2d_z
!deck print_psi_3d_z
!***begin prologue     print_psi_3d_z
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***references
!***routines called
!***end prologue       print_psi_3d_z
!
  SUBROUTINE print_psi_3d_z(wave_function)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                        :: i, j, k
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: wave_function
!
  write(iout,1) title
  write(iout,2)
  DO i = 1,nphy(1)
     DO j = 1, nphy(2)
        DO k = 1, nphy(3)
           write(iout,3) grid(1)%pt(i), grid(2)%pt(j), grid(3)%pt(k),  &
                         wave_function(k,j,i,1), wave_function(k,j,i,2)
        END DO
     END DO
  END DO
1 FORMAT(a80)
2 FORMAT('     x      ',5x,'     y      ',5x,'     z     ','         &
                                         Real Psi  ',5x,' Imag Psi ')
3 FORMAT(e12.6,5x,e12.6,5x,e12.6,5x,e12.6,5x,e12.6)
END SUBROUTINE print_psi_3d_z
END MODULE plot_wavefunction

