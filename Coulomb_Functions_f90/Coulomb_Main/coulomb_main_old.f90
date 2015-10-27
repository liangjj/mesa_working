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
  USE Coulomb_Module
  USE Auxilliary_functions_Module
  USE Series_Expansion
  USE Asymptotic_Expansion
  IMPLICIT NONE
  LOGICAL                                  :: dollar
  REAL*8                                   :: fpkey
  INTEGER                                  :: input
  INTEGER                                  :: output
  INTEGER                                  :: intkey
  INTEGER                                  :: i
  INTEGER                                  :: j
  LOGICAL                                  :: logkey
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
     quantities_returned = chrkey(card,'quantities_returned',          &
                                  'functions_and_derivatives',' ')
     type = chrkey(card,'function_type','coulomb',' ')
     r_series = fpkey(card,'series_cut_off',10.d0,' ')
     r_asymptotic = fpkey(card,'asymptotic_cut_off',10.d0,' ')
     print_sigma_l = logkey(card,'print_sigma_l',.false.,' ')
     print_long_range_coefficients = logkey(card,'print_long_range_coefficients',    &
                                            .false.,' ')
     print_convergence = logkey(card,'print_convergence',.false.,' ')
     print_short_range_coefficients = logkey(card,'print_short_range_coefficients',  &
                                             .false.,' ')
     series_size = intkey(card,'number_of_terms_in_power_series',50,' ')
     asymptotic_size = intkey(card,'number_of_terms_in_asymptotic_series',50,' ')
     ALLOCATE(rho(1:number_of_r_values), rho_inv(1:number_of_r_values))
     ALLOCATE(fc(1:l_max+1), dfc(1:l_max+1), gc(1:l_max+1), dgc(1:l_max+1))
     rho(1) = 0.d0
     DO i=2, number_of_r_values
        rho(i) = rho(i-1) + r_step
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
     write(iout,4) print_sigma_l, print_long_range_coefficients,       &
                   print_convergence, print_short_range_coefficients
!     CALL short_range_coefficients
!     CALL long_range_coefficients
     IF(energy_class == 'positive' ) THEN  
        write(iout,5)
        DO i = 1, number_of_r_values
           r = rho(i)
           r_inv = rho_inv(i)
!           IF ( r <= r_series) THEN
!                WRITE(iout,6)
!**********************************************************************c
!            use series expansion for small rho                        c
!**********************************************************************c
!                CALL series_expansion_regular_positive_energy_function
!                CALL series_expansion_irregular_positive_energy_function
           Call Coulfg (rho(i),eta_in,fc,gc,dfc,dgc)
           IF (ifail /= 0) THEN
               write(iout,*) 'Warning message !!!'
           END IF
           write(iout,4) rho(i)
           DO j = l_min+1 , l_max+1
              fl  = fc(j)
              dfl = dfc(j)
              gl  = gc(j)
              dgl = dgc(j)
              wronskian = fl*dgl - dfl*gl
              write(iout,5) j-1, fl, dfl , gl , dgl, wronskian
           END DO
!           ELSE
!                WRITE(iout,8)
!                CALL asymptotic_expansion_positive_energy_function
!                wronskian = fl*dgl-dfl*gl 
!                WRITE(iout,7) r, fl, dfl, gl, dgl, wronskian
!           END IF
        END DO
!     ELSE IF(energy_class == 'negative') THEN
!        IF ( r <= r_series) THEN
!            write(iout,5)
!        ELSE IF( r >= r_asymptotic) THEN
!            write(iout,10)
!        END IF
!        DO i = 1, number_of_r_values
           r = rho(i)
           r_inv = rho_inv(i)
           IF ( r <= r_series) THEN
                write(iout,6)
!**********************************************************************c
!             use series expansion for small rho                        c
!**********************************************************************c
                CALL series_expansion_regular_negative_energy_function
                CALL series_expansion_irregular_negative_energy_function
                wronskian = fl*dgl-dfl*gl 
                WRITE(iout,7) r, fl, dfl, gl, dgl, wronskian
           ELSE
                write(iout,8)
                CALL asymptotic_expansion_irregular_negative_energy_function
                WRITE(iout,7) r, gl, dgl
           END IF
        END DO
     END IF
  END IF
  WRITE(iout,1)
1 FORMAT('************************************************************************')
2 FORMAT(25X,'Coulomb Functions')
3 FORMAT(/,5X,'Energy                    = ',e15.8,1X,'k(kappa)      = ',e15.\
8,  &
         /,5x,'Smallest Angular Momentum = ',i3,                             \
    &
           1x,'Largest Angular Momentum  = ',i3,/,5x,'Eta  = ',e15.8)
4 FORMAT(/,5x,'print_sigma_l = ',l1,1x,'print_long_range_coefficients = ',l1   &
         /,5x,'print_convergence = 'l1,1x,'print_short_range_coefficients = ',l1 )
5 FORMAT(/,11x,'r', 13x,'R_l',13x,'DR_l',13x,'G_l',12x,'DG_l',10x,'Wron',      &
         /,11x, '-',13x,'---',13x,'----',13x,'---',12x,'----',10x,'----')
6 FORMAT(35x,'Series Expansion')
7 FORMAT(3x,6(1x,e15.8))
8 FORMAT(35x,'Asymptotic Expansion')
9 FORMAT (/,5x,'Radial Point out of Range')
10 FORMAT(/,11x,'r', 13x,'G_l',13x,'DG_l')
END PROGRAM Coulomb_Main
