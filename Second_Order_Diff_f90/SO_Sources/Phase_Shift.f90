!deck Phase_Shift
!***begin prologue     Phase_Shift
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           space, propagation, dvr, finite difference
!***author             schneider, b. i.(nsf)
!***source
!***purpose            
!***                   
!***                   
!***                   
!***                   
!***                   
!***                   
!***                   
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
  Subroutine Phase_Shift (surface_vector,eigenvalues,R_Box,energy,n)
  USE input_output
  USE dvr_shared,              ONLY : charge  
  USE Coulomb_Functions_Module,                                                    &
                               ONLY : eye, one, k, angular_momentum, eta_in, fgd,  &
                                      c_arg, sigma_l, cgamma, type,                &
                                      quantities_returned, pi, f, g,               &
                                      f_d, g_d, asymptotic_size,                   &
                                      asymptotic_series, print_sigma_l,            &
                                      print_long_range_coefficients,               &
                                      print_convergence,                           &
                                      fl, dfl, gl, dgl, wronskian,                 &
                                      phase_factor, r, r_inv,                      &
                                      s_sigma, c_sigma,                            &
                                      asymptotic_expansion_positive_energy_function 
  USE Asymptotic_Module,       ONLY : positive_energy_long_range_coefficients
                                             
  IMPLICIT NONE
  INTEGER                                  :: n
  REAL*8, DIMENSION(n)                     :: surface_vector
  REAL*8, DIMENSION(n)                     :: eigenvalues
  REAL*8                                   :: R_mat
  REAL*8                                   :: R_Box
  REAL*8                                   :: j
  REAL*8                                   :: dj
  REAL*8                                   :: ddj
  REAL*8                                   :: y
  REAL*8                                   :: dy
  REAL*8                                   :: ddy
  REAL*8                                   :: s_n
  REAL*8                                   :: c_n
  REAL*8                                   :: d_r
  REAL*8                                   :: tan_delta
  REAL*8                                   :: phase
  REAL*8                                   :: energy
  REAL*8                                   :: ratio
  REAL*8                                   :: s_coul
  REAL*8                                   :: c_coul
  INTEGER                                  :: ene
!
  k=sqrt(2.d0*energy)
  R_mat = 0.0d0
  DO ene=1,n
     R_mat = R_mat + surface_vector(ene) * surface_vector(ene)   &
                                          /                      &
                             ( energy - eigenvalues(ene) )
  END DO
  R_mat = - .5d0 * R_mat
  r = k*R_Box
  r_inv = 1.d0/r
  IF (charge == 0.d0) THEN
!
!      j goes like sin(kr - l*pi/2), y like  -cos(kr-l*pi/2)
!
      Call rbes('ricatti-bessel',angular_momentum,r,j,dj,ddj,y,dy,ddy)
      wronskian = j*dy - y*dj      
      write(iout,1) j, dj, y, dy, wronskian
!
!     Change to cos not minus cos to make consistent with coulomb
!
      y = - y
      dy = - dy
      dj = k * dj
      dy = k * dy
      ratio = j /dj
  ELSE
      eta_in = charge/k
      ALLOCATE(asymptotic_series(angular_momentum)%a_i(0 : asymptotic_size),       &
               asymptotic_series(angular_momentum)%b_i(0 : asymptotic_size) )
      Call positive_energy_long_range_coefficients(                                &
                   asymptotic_series(angular_momentum)%a_i(0 : asymptotic_size),   &
                   asymptotic_series(angular_momentum)%b_i(0 : asymptotic_size) )
      sigma_l=c_arg(cgamma(angular_momentum + one + eye*eta_in))
      type='coulomb'
      quantities_returned='functions_and_derivatives'
      Call asymptotic_expansion_positive_energy_function                           &
                                             (asymptotic_series(angular_momentum)%a_i,    &
                                              asymptotic_series(angular_momentum)%b_i )
      wronskian = fl*dgl - gl*dfl
      write(iout,1) fl, dfl, gl, dgl, wronskian
!
!     We now write the wavefunction as a linear combination of cos(sigma_l)
!     and sin(sigma_l) and compute sigma_l
!
      j = phase_factor(1)
      dj = k*phase_factor(3)
      y = phase_factor(2)
      dy = k*phase_factor(4)
      DEALLOCATE(asymptotic_series(angular_momentum)%a_i, asymptotic_series(angular_momentum)%b_i )
  END IF
  tan_delta = ( (R_mat * dj - j) / (y - R_mat * dy) )
  phase = atan(tan_delta)
  write(iout,2) charge, energy, R_mat, sigma_l, phase
1 Format(/,10x,'Regular Function   = ',e15.8,2x,            &
               'Derivative of Regular Function   = ',e15.8, &
         /,10x,'Irregular Function = ',e15.8,2x,            &
               'Derivative of Irregular Function = ',e15.8, &
         /,10x,'Wronskian = ',e15.8)
2 Format(/,10x,'Charge   = ',e15.8,                         &
         /,10x,'Energy   = ',e15.8,                         &
         /,10x,'R Matrix = ',e15.8,                         &
         /,10x,'Analytic Value of Phase Shift   = ',e15.8,  &
         /,10x,'Calculated Value of Phase Shift = ',e15.8)
END Subroutine Phase_Shift
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
