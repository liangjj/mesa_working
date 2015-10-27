!deck @(#)Coulomb_Main
!***begin prologue     Coulomb_Main
!***date written       920525   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           m7000, link 7000, spline
!***author             schneider, b. (nsf)
!***source             m7004
!***purpose            regular and irregular coulomb functions
!***
!***description        regular and irregular coulomb function are computed
!***                   using a series expansion at small r and an asymptotic
!***                   expansion at large r.
!***references

!***routines called    iosys, util and mdutil
!***end prologue       Coulomb_Main

  Program Coulomb_Main
  USE Coulomb_Variables_Module
  USE Special_functions_Module
  USE Series_Module
  USE Asymptotic_Module
  USE Coulomb_Functions_Module
  IMPLICIT NONE
  LOGICAL                                  :: dollar
  REAL*8                                   :: fpkey
  REAL*8                                   :: mult_fact
  INTEGER                                  :: input
  INTEGER                                  :: output
  INTEGER                                  :: intkey
  INTEGER                                  :: pt
  LOGICAL                                  :: logkey
  CHARACTER (LEN=80)                       :: chrkey
  CHARACTER (LEN=16)                       :: fptoc
  COMMON /io/ input, output
!
  input = inp
  output = iout

!
!  Open the input and output files
!
  OPEN(input,file='Coulomb.inp',status='old')
  OPEN(output,file='Coulomb.out',status='unknown')
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
  ALLOCATE (fact(0:max_fact))
  CALL factorials
  IF ( dollar('$coulomb',card,cpass,inp) ) then
!
!**********************************************************************c
!                calculate k or kappa and the value of the             c
!                             coulomb eta                              c
!**********************************************************************c
     charge = fpkey(card,'charge',-1.d0,' ')
     energy = fpkey(card,'energy',1.d0,' ')
     l_min = intkey(card,'smallest_angular_momentum',0,' ')
     l_max = intkey(card,'largest_angular_momentum',0,' ')
     r_min = fpkey(card,'smallest_r',0.d0,' ')
     r_max = fpkey(card,'largest_r',0.d0,' ')
     r_step = fpkey(card,'delta_r',0.d0,' ')
     number_of_r_values = ( r_max - r_min ) / r_step + 1
     quantities_returned = chrkey(card,'quantities_returned',                  &
                                  'functions_and_derivatives',' ')
     type = chrkey(card,'function_type','coulomb',' ')
     r_series = fpkey(card,'series_cut_off',.05d0,' ')
     r_asymptotic = fpkey(card,'asymptotic_cut_off',15.d0,' ')
     print_sigma_l = logkey(card,'print_sigma_l',.false.,' ')
     print_long_range_coefficients = logkey(card,                              &
                                           'print_long_range_coefficients',    &
                                            .false.,' ')
     print_convergence = logkey(card,'print_convergence',.false.,' ')
     print_short_range_coefficients = logkey(card,                             &
                                             'print_short_range_coefficients', &
                                             .false.,' ')
     series_size = intkey(card,'number_of_terms_in_power_series',50,' ')
     asymptotic_size = intkey(card,'number_of_terms_in_asymptotic_series',50,' ')
  ELSE
     write(iout,*) 'STOP INPUT CARD ERROR'
     stop
  END IF
  ALLOCATE(rho(1:number_of_r_values), rho_inv(1:number_of_r_values))
  ALLOCATE(fgd(1:l_max+1,9))
  rho(1) = r_min
  DO pt=2, number_of_r_values
     rho(pt) = rho(pt-1) + r_step
  END DO
  IF(energy >= zero ) THEN
     energy_class = 'positive'
     k = SQRT(energy*two)
  ELSE
     energy_class = 'negative'
     k = SQRT ( two * abs(energy) )
  END IF
  eta_in=charge/k
  rho(:) = k * rho(:)
  rho_inv(1:number_of_r_values) = 1.d0 / rho(1:number_of_r_values)
  WRITE(iout,3) energy, k, l_min, l_max, eta_in
!**********************************************************************c
!        calculate the coefficients for the series expansion           c
!        at small rho and the asymptotic expansion at large rho        c
!**********************************************************************c
  ALLOCATE( power_series(0:l_max), asymptotic_series(0:l_max) )
  IF ( energy_class == 'positive' ) THEN
       DO angular_momentum = l_min, l_max
          ALLOCATE(                                                                                          &
          power_series(angular_momentum)%a(angular_momentum + 1 : angular_momentum + 1 + series_size),       &
          power_series(angular_momentum)%b(-angular_momentum - 1 : - angular_momentum  + series_size) )                    
          Call positive_energy_short_range_coefficients(power_series(angular_momentum)%a,                    &
                                                        power_series(angular_momentum)%b )
          ALLOCATE(                                                                                          &
                   asymptotic_series(angular_momentum)%a_i(0 : asymptotic_size),                             &
                   asymptotic_series(angular_momentum)%b_i(0 : asymptotic_size) )
          Call positive_energy_long_range_coefficients(                                                      &
                                       asymptotic_series(angular_momentum)%a_i(0 : asymptotic_size),         &
                                       asymptotic_series(angular_momentum)%b_i(0 : asymptotic_size) )  
       END DO
  ELSE IF ( energy_class == 'negative' ) THEN
       DO angular_momentum = l_min, l_max
          ALLOCATE(                                                                                          &
          power_series(angular_momentum)%a_0(0 : series_size),                                               &
          power_series(angular_momentum)%b_0(0 : series_size),                                               &
          power_series(angular_momentum)%c_0(0 : series_size),                                               &
          power_series(angular_momentum)%d_0(0 : series_size)  )
          Call negative_energy_short_range_coefficients(power_series(angular_momentum)%a_0,                  &
                                                        power_series(angular_momentum)%b_0,                  &
                                                        power_series(angular_momentum)%c_0,                  &
                                                        power_series(angular_momentum)%d_0 )
          ALLOCATE( asymptotic_series(angular_momentum)%e_0(0 : asymptotic_size)  )
          Call negative_energy_long_range_coefficients(                                                      &
                                       asymptotic_series(angular_momentum)%e_0(0 : asymptotic_size) )  
       END DO
  END IF
  write(iout,4) print_sigma_l, print_long_range_coefficients,                                                &
                print_convergence, print_short_range_coefficients
  IF ( energy_class == 'positive' ) THEN
       DO pt = 1, number_of_r_values
          r = rho(pt)
          r_inv = rho_inv(pt)
          IF (r <= r_series) THEN
              Write(iout,8)
               DO angular_momentum=l_min,l_max
                  Call series_expansion_regular_positive_energy_function(power_series(angular_momentum)%a)
                  fgd(angular_momentum+1,1) = fl
                  fgd(angular_momentum+1,3) = dfl
                  Call series_expansion_irregular_positive_energy_function(power_series(angular_momentum)%a, &
                                                                           power_series(angular_momentum)%b, &
                                                                           .true. )
                  fgd(angular_momentum+1,2) = gl
                  fgd(angular_momentum+1,4) = dgl
              END DO
          ELSE IF (r > r_series.and. r<r_asymptotic) THEN
                  Write(iout,9)
                  Call Coulfg (rho(pt),eta_in,fgd(:,1),fgd(:,2),fgd(:,3),fgd(:,4))
                  IF (ifail /= 0) THEN
                      write(iout,*) 'Warning message !!!'
                  END IF
                  IF(iexp.gt.int_one) THEN
                     mult_fact = 10**iexp
                     fgd(:,1) = fgd(:,1) / mult_fact
                     fgd(:,3) = fgd(:,3) / mult_fact
                     fgd(:,2) = fgd(:,2) * mult_fact
                     fgd(:,4) = fgd(:,4) * mult_fact
                  END IF
          ELSE IF (r >= r_asymptotic) THEN
               Write(iout,10)
               DO angular_momentum=l_min,l_max
                  Call asymptotic_expansion_positive_energy_function                                          &
                                              (asymptotic_series(angular_momentum)%a_i,                       &
                                               asymptotic_series(angular_momentum)%b_i )
                  fgd(angular_momentum+1,1) = fl
                  fgd(angular_momentum+1,3) = dfl
                  fgd(angular_momentum+1,2) = gl
                  fgd(angular_momentum+1,4) = dgl
                  fgd(angular_momentum+1,6) = f
                  fgd(angular_momentum+1,7) = g
                  fgd(angular_momentum+1,8) = f_d
                  fgd(angular_momentum+1,9) = g_d
               END DO
          END IF
          title='Coulomb Functions at rho = '//fptoc(rho(pt))
          fgd(:,5) = fgd(:,1)*fgd(:,4) - fgd(:,2)*fgd(:,3)
          write(iout,5) title
          write(iout,6)
          DO angular_momentum=l_min, l_max
             write(iout,7) angular_momentum, fgd(angular_momentum+1,:)
          END DO
      END DO
  ELSE IF ( energy_class == 'negative' ) THEN
       DO pt = 1, number_of_r_values
          r = rho(pt)
          r_inv = rho_inv(pt)
          IF (r <= r_series) THEN
              Write(iout,8)
              DO angular_momentum=l_min,l_max
                 Call series_expansion_regular_negative_energy_function                                       &
                                              (power_series(angular_momentum)%a_0)
                 fgd(angular_momentum+1,1) = fl
                 fgd(angular_momentum+1,3) = dfl
                 Call series_expansion_irregular_negative_energy_function(                                    &
                                                     power_series(angular_momentum)%a_0,                      &
                                                     power_series(angular_momentum)%b_0,                      &
                                                     power_series(angular_momentum)%c_0 )
                 fgd(angular_momentum+1,2) = gl
                 fgd(angular_momentum+1,4) = dgl
              END DO
          ELSE IF (r > r_series.and. r<r_asymptotic) THEN
                   Write(iout,9)
                   Call Coulfg (rho(pt),eta_in,fgd(:,1),fgd(:,2),fgd(:,3),fgd(:,4))
                   IF (ifail /= 0) THEN
                       write(iout,*) 'Warning message !!!'
                   END IF
          ELSE IF (r >= r_asymptotic) THEN
               Write(iout,10)
               DO angular_momentum=l_min,l_max
                  Call asymptotic_expansion_irregular_negative_energy_function                                &
                                              (asymptotic_series(angular_momentum)%e_0  )
                  fgd(angular_momentum+1,1) = zero
                  fgd(angular_momentum+1,3) = zero                
                  fgd(angular_momentum+1,2) = gl
                  fgd(angular_momentum+1,4) = dgl
               END DO
          END IF
          title='Coulomb Functions at rho = '//fptoc(rho(pt))
          fgd(:,5) = fgd(:,1)*fgd(:,4) - fgd(:,2)*fgd(:,3)
          write(iout,5) title
          write(iout,6)
          DO angular_momentum=l_min, l_max
             write(iout,7) angular_momentum, fgd(angular_momentum+1,:)
          END DO
       END DO
  END IF

  WRITE(iout,1)
1 FORMAT('************************************************************************')
2 FORMAT(25X,'Coulomb Functions')
3 FORMAT(/,5X,'Energy                    = ',e15.8,1X,'k(kappa)      = ',e15.8,  &
         /,5x,'Smallest Angular Momentum = ',i3,                                 &
           1x,'Largest Angular Momentum  = ',i3,/,5x,'Eta  = ',e15.8)
4 FORMAT(/,5x,'print_sigma_l = ',l1,1x,'print_long_range_coefficients = ',l1     &
         /,5x,'print_convergence = 'l1,1x,'print_short_range_coefficients = ',l1 )
5 FORMAT(/,30x,a80)
6 FORMAT(/,7x,'L',8x,'F',14x,'G',13x,'DF',12x,'DG',13x,'WRON')
7 FORMAT(/,5x,i3,5e15.8)
8 Format(/,5x,'Series Expansion')
9 Format(/,5x,'Steeds Method')
10 Format(/,5x,'Asymptotic Expansion')
END PROGRAM Coulomb_Main
