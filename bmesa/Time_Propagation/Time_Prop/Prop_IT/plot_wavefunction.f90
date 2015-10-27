!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        MODULE plot_wavefunction
                        USE dvrprop_global
                        USE dvr_shared  
                        USE dvr_global
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
             MODULE PROCEDURE plot_psi_d,                               &
                              plot_psi_z
                    END INTERFACE plot_psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        INTERFACE print_psi
             MODULE PROCEDURE print_psi_d,                              &
                              print_psi_z
                    END INTERFACE print_psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck plot_psi_d.f
!***begin prologue     plot_psi_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_psi_d
  SUBROUTINE plot_psi_d(vector,t)
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)                 :: vector
  REAL*8                                 :: rtemp
  INTEGER                                :: t, i, j, k, iostat, count
  CHARACTER (LEN=16)                     :: fptoc
  IF (log_main(8)) THEN
      title='solution at t = '//fptoc(tim_pts(t+1))
      CALL print_psi(vector)
  END IF
  IF (plot) THEN
      OPEN (UNIT=iplot(3),FILE='wavefunction-'//fptoc(tim_pts(t+1)),         &
            ACCESS='sequential',FORM='formatted',                            &
            IOSTAT=IOSTAT,STATUS='unknown')
      IF(IOSTAT /= 0) THEN
         CALL lnkerr('error in file handling')
      END IF
      OPEN (UNIT=iplot(4),FILE='modulus_wavefunction-'//fptoc(tim_pts(t+1)),  &
            ACCESS='sequential',FORM='formatted',                             &
            IOSTAT=IOSTAT,STATUS='unknown')
      IF(IOSTAT /= 0) THEN
         CALL lnkerr('error in file handling')
      END IF
      IF (spdim == 1 ) THEN
          DO i = 1,nphy(1)
             rtemp = ABS ( vector(i)  )
             write(iplot(3),*) grid(1)%pt(i), vector(i)
             write(iplot(4),*) grid(1)%pt(i), rtemp
          END DO
      ELSE IF (spdim == 2 ) THEN
          count = 0
          DO i = 1,nphy(2)
             DO j = 1,nphy(1)
                count = count + 1
                rtemp = ABS ( vector(count) )
                write(iplot(3),*) grid(2)%pt(i), grid(1)%pt(j), vector(count)
                write(iplot(4),*) grid(2)%pt(i), grid(1)%pt(j), rtemp
             END DO
          END DO
      ELSE IF ( spdim == 3 ) THEN
          count = 0
          DO i = 1,nphy(3)
             DO j = 1,nphy(2)
                DO k = 1,nphy(1)
                   count = count + 1
                   rtemp = ABS ( vector(count) )
                   write(iplot(3),*) grid(3)%pt(i), grid(2)%pt(j),         &
                                     grid(1)%pt(k), vector(count)      
                   write(iplot(4),*) grid(3)%pt(i), grid(2)%pt(j),         &
                                     grid(3)%pt(k), rtemp
                END DO
             END DO
          END DO
      END IF
  END IF
END SUBROUTINE plot_psi_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck plot_psi_z.f
!***begin prologue     plot_psi_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_psi_z
  SUBROUTINE plot_psi_z(vector,t)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d)             :: vector
  REAL*8                                 :: rtemp
  INTEGER                                :: t, i, j, k, iostat, count
  CHARACTER (LEN=16)                     :: fptoc
  IF (log_main(8)) THEN
      title='solution at t = '//fptoc(tim_pts(t+1))
      CALL print_psi(vector)
  END IF
  IF (plot) THEN
      OPEN (UNIT=iplot(3),FILE='wavefunction-'//fptoc(tim_pts(t+1)),         &
            ACCESS='sequential',FORM='formatted',                            &
            IOSTAT=IOSTAT,STATUS='unknown')
      IF(IOSTAT /= 0) THEN
         CALL lnkerr('error in file handling')
      END IF
      OPEN (UNIT=iplot(4),FILE='modulus_wavefunction-'//fptoc(tim_pts(t+1)),  &
            ACCESS='sequential',FORM='formatted',                             &
            IOSTAT=IOSTAT,STATUS='unknown')
      IF(IOSTAT /= 0) THEN
         CALL lnkerr('error in file handling')
      END IF
      IF (spdim == 1 ) THEN
          DO i = 1,nphy(1)
             rtemp = SQRT ( vector(i) * conjg( vector(i) ) )
             write(iplot(3),*) grid(1)%pt(i), vector(i)
             write(iplot(4),*) grid(1)%pt(i), rtemp
          END DO
      ELSE IF (spdim == 2 ) THEN
          count = 0
          DO i = 1,nphy(2)
             DO j = 1,nphy(1)
                count = count + 1
                rtemp = SQRT ( vector(count) * conjg( vector(count) ) )
                write(iplot(3),*) grid(2)%pt(i), grid(1)%pt(j), vector(count)
                write(iplot(4),*) grid(2)%pt(i), grid(1)%pt(j), rtemp
             END DO
          END DO
      ELSE IF ( spdim == 3 ) THEN
          count = 0
          DO i = 1,nphy(3)
             DO j = 1,nphy(2)
                DO k = 1,nphy(1)
                   count = count + 1
                   rtemp = SQRT ( vector(count) * conjg( vector(count) ) )
                   write(iplot(3),*) grid(3)%pt(i), grid(2)%pt(j),         &
                                     grid(1)%pt(k), vector(count)      
                   write(iplot(4),*) grid(3)%pt(i), grid(2)%pt(j),         &
                                     grid(3)%pt(k), rtemp
                END DO
             END DO
          END DO
      END IF
  END IF
END SUBROUTINE plot_psi_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck print_psi_d
!***begin prologue     print_psi_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***references
!***routines called
!***end prologue       print_psi_d
!
  SUBROUTINE print_psi_d(vector)
  IMPLICIT NONE
  INTEGER                                :: i, j, k, count
  REAL*8, DIMENSION(n3d)                 :: vector
!
  write(iout,1) title
  IF ( spdim == 1 ) THEN
     write(iout,2)
     DO i=1,nphy(1)
        write(iout,5) grid(1)%pt(i), vector(i)
     END DO
  ELSE IF ( spdim == 2 ) THEN
     write(iout,3)
     count = 0 
     DO i = 1,nphy(1)
        DO j = 1, nphy(2)
           count = count + 1 
           write(iout,6) grid(1)%pt(i), grid(2)%pt(j),                 &
                         vector(count)
        END DO
     END DO
  ELSE IF ( spdim == 3 ) THEN
     write(iout,4)
     count = 0 
     DO i = 1,nphy(1)
        DO j = 1, nphy(2)
           DO k = 1, nphy(3)
              count = count + 1 
              write(iout,7) grid(1)%pt(i), grid(2)%pt(j), grid(3)%pt(k),  &
                            vector(count)
           END DO
        END DO
     END DO
  END IF
1 FORMAT(a80)
2 FORMAT('     x      ',5x,'   Psi  ')
3 FORMAT('     x      ',5x,'     y      ',5x,'  Psi  ')
4 FORMAT('     x      ',5x,'     y      ',5x,'     z     ','         &
                                         Psi  ')
5 FORMAT(e12.6,5x,e12.6)
6 FORMAT(e12.6,5x,e12.6,5x,e12.6)
7 FORMAT(e12.6,5x,e12.6,5x,e12.6,5x,e12.6)
END SUBROUTINE print_psi_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE print_psi_z(vector)
  IMPLICIT NONE
  INTEGER                                :: i, j, k, count
  COMPLEX*16, DIMENSION(nphy(1))         :: vector
!
  write(iout,1) title
  IF ( spdim == 1 ) THEN
     write(iout,2)
     DO i=1,nphy(1)
        write(iout,5) grid(1)%pt(i), vector(i)
     END DO
  ELSE IF ( spdim == 2 ) THEN
     write(iout,3)
     count = 0 
     DO i = 1,nphy(1)
        DO j = 1, nphy(2)
           count = count + 1 
           write(iout,6) grid(1)%pt(i), grid(2)%pt(j),                 &
                         vector(count)
        END DO
     END DO
  ELSE IF ( spdim == 3 ) THEN
     write(iout,4)
     count = 0 
     DO i = 1,nphy(1)
        DO j = 1, nphy(2)
           DO k = 1, nphy(3)
              count = count + 1 
              write(iout,7) grid(1)%pt(i), grid(2)%pt(j), grid(3)%pt(k),  &
                            vector(count)
           END DO
        END DO
     END DO
  END IF
1 FORMAT(a80)
2 FORMAT('     x      ',5x,'   Real Psi  ',5x,'Imag Psi ')
3 FORMAT('     x      ',5x,'     y      ',5x,'  Real Psi  ',5x,' Imag Psi ')
4 FORMAT('     x      ',5x,'     y      ',5x,'     z     ','         &
                                         Real Psi  ',5x,' Imag Psi ')
5 FORMAT(e12.6,5x,e12.6,5x,e12.6)
6 FORMAT(e12.6,5x,e12.6,5x,e12.6,5x,e12.6)
7 FORMAT(e12.6,5x,e12.6,5x,e12.6,5x,e12.6,5x,e12.6)
END SUBROUTINE print_psi_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE plot_wavefunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

