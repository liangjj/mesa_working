!*********************************************************************
                      MODULE Special_Functions_Module
                      USE Coulomb_Variables_Module
                      IMPLICIT NONE
  REAL*8                                             :: eta_sq
  REAL*8                                             :: two_pi_eta
  REAL*8                                             :: sum_1
  REAL*8                                             :: sum_2
  REAL*8                                             :: denominator
  REAL*8, DIMENSION(:),   ALLOCATABLE                :: fact
  COMPLEX*16                                         :: eta_fac
  COMPLEX*16                                         :: com_eta
  COMPLEX*16                                         :: arg
  COMPLEX*16                                         :: num
  COMPLEX*16                                         :: cr_l
  INTEGER                                            :: i
  INTEGER                                            :: max_fact = 100
  INTEGER                                            :: two_el
  INTEGER                                            :: two_el_1
  INTEGER                                            :: two_el_m_1
  INTEGER                                            :: l_1
  INTEGER                                            :: two_i
  INTEGER                                            :: two_i_1
  INTEGER                                            :: two_i_m_1
  INTEGER                                            :: n_terms
  INTEGER                                            :: l_fac
  INTEGER                                            :: i_test
!=======================================================================
!=======================================================================
                        Contains
!=======================================================================
!=======================================================================
!deck factorials
!***begin prologue     factorials
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!                      
!***description 
!***references 
!***routines called
!***end prologue       factorials
  Subroutine factorials
  IMPLICIT NONE
  fact(0) = one
  DO i = 1, max_fact
     fact(i) = i * fact(i-1)
  END DO
  END SUBROUTINE factorials
!=======================================================================
!=======================================================================
!deck c0_sq_fun
!***begin prologue     c0_sq_fun
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!                                         2
!***description        coulomb function c
!                                        0
!***references         NBS handbook

!***routines called

!***end prologue       c0_sq_fun
  REAL*8 FUNCTION c0_sq_fun()
  IMPLICIT NONE
  two_pi_eta = two * pi * eta_in
  c0_sq_fun = two_pi_eta/(EXP(two_pi_eta) - one)
  END FUNCTION c0_sq_fun
!=======================================================================
!=======================================================================
!deck  cl_fun
!***begin prologue     cl_fun
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            compute c  function for positive energy coulomb
!***                            l                           functions.
!***description
!***references         NBS handbook
!***routines called    cgamma(complex gamma function:clams)
!***end prologue       m2705
  REAL*8 Function cl_fun()
  IMPLICIT NONE
  cl_fun = (two**angular_momentum)*EXP(-pi*eta_in*half)
  l_1 = angular_momentum + 1
  two_el = angular_momentum + angular_momentum
  arg = l_1 + eye*eta_in
  cl_fun=cl_fun*ABS(cgamma(arg))/fact( two_el + 1 )
  END FUNCTION cl_fun
!=======================================================================
!=======================================================================
!deck dl_fun
  REAL*8 FUNCTION dl_fun()
!***begin prologue     dl_fun
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***description        coulomb d  function
!                               l
!                      this function is trivially related to c
!                                                             l
!***references         NBS handbook

!***routines called    cl_fun(code)
!***end prologue       dl_fun
  IMPLICIT NONE
  dl_fun = one /( (angular_momentum + angular_momentum + 1)*cl_fun() )
  END FUNCTION dl_fun
!=======================================================================
!=======================================================================
!deck pl_fun
  REAL*8 FUNCTION pl_fun()
!***begin prologue     pl_fun
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***description        coulomb p  function
!                               l
!***references         NBS handbook
!***routines called
!***end prologue       pl_fun
  IMPLICIT NONE
  pl_fun = two*eta_in
  IF (angular_momentum > 0) THEN
      eta_sq = eta_in*eta_in
      pl_fun = pl_fun*(one+eta_sq)/three
      IF (angular_momentum > 1) THEN
          DO i = 2, angular_momentum
             two_i = i + i
             two_i_1 = two_i + 1
             two_i_m_1 = two_i -1
             pl_fun = pl_fun*four*(angular_momentum*angular_momentum + eta_sq)       &
                                /                              &
                      (two_i_1*two_i*two_i*two_i_m_1)
          END DO
      END IF
  END IF
  END FUNCTION pl_fun
!=======================================================================
!=======================================================================
!deck ql_pl_fun
!***begin prologue     ql_pl_fun
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***description        q /p  function for construction of irregular positive
!                       l  l energy coulomb functions
!***references          nbs mathematical handbook (abromowitz and stegun)
!***routines called    cpsi(complex psi function:clams)
!***end prologue       ql_pl_fun
  REAL*8 FUNCTION ql_pl_fun()
  IMPLICIT NONE
  eta_fac = one + eta_in*eye
  IF (angular_momentum == 0) THEN
      ql_pl_fun = - one + REAL( cpsi(eta_fac) ) + two*eulerc
  ELSE
      eta_sq = eta_in*eta_in
      sum_1 = zero
      DO i = 1 , angular_momentum
         sum_1 = sum_1 + i/( i * i + eta_sq)
      END DO
      sum_2 = zero
      DO i = 1 ,angular_momentum + angular_momentum + 1
         sum_2 = sum_2 + one / i
      END DO
      ql_pl_fun = sum_1 - sum_2 + REAL( cpsi(eta_fac) ) + two*eulerc +  &
                               rl_fun()/pl_fun()
  END IF
  END FUNCTION ql_pl_fun
!=======================================================================
!=======================================================================
!deck rl_fun
!***begin prologue     rl_fun
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***description        r function for construction of irregular positive
!                       l energy coulomb functions
!***references          nbs mathematical handbook (abromowitz and stegun)
!***routines called
!***end prologue       rl_fun
  REAL*8 FUNCTION rl_fun()
  IMPLICIT NONE
  IF (angular_momentum == 0) THEN
      rl_fun = zero
  ELSE
      n_terms = angular_momentum+angular_momentum+1
      l_fac = n_terms
      com_eta = eta_in*eye - angular_momentum
      num = one
      cr_l = num/l_fac
      DO i= 2 , n_terms
         l_fac = l_fac-1
         num = num * two * com_eta
         denominator = l_fac*fact(i-1)
         cr_l = cr_l + num / denominator
         com_eta = com_eta + one
      END DO
      rl_fun = imag(cr_l)
      i_test = angular_momentum - 2 * (angular_momentum/2)
      IF (i_test == 0) THEN
          rl_fun = -rl_fun
      END IF
      rl_fun = rl_fun / fact(angular_momentum + angular_momentum)
  END IF
  END FUNCTION rl_fun
!====================================================================
!====================================================================
!deck a_fun
!***begin prologue     a_fun
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***description        series expansion coefficients for regular positive energy coulomb functions
!                      
!***references          nbs mathematical handbook (abromowitz and stegun)
!***routines called
!***end prologue       a_fun
  Subroutine a_fun(a)
  IMPLICIT NONE
  INTEGER                 :: upper 
  REAL*8, DIMENSION(:)    :: a(angular_momentum + 1 : angular_momentum + 1 + series_size)
  REAL*8                  :: argument
  INTEGER                 :: i 
  INTEGER                 :: l_one 
  INTEGER                 :: l_two 
  INTEGER                 :: l_three 
  INTEGER                 :: l_upper 
  l_one = angular_momentum + int_one
  l_two = l_one + int_one
  l_three = l_two + int_one
  l_upper = l_one + series_size
  argument = two * eta_in
  a(l_one) = one
  a(l_two)= eta_in / l_one
  DO i = l_three , l_upper 
     a(i) = argument * a(i-1) - a(i-2)
     a(i) = a(i)/ ( (i + angular_momentum ) * ( i - l_one) )
  END DO
  END Subroutine a_fun
!====================================================================
!====================================================================
!deck b_fun
!***begin prologue     b_fun
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***description        series expansion coefficients for irregular positive energy coulomb functions
!                      
!***references          nbs mathematical handbook (abromowitz and stegun)
!***routines called
!***end prologue       b_fun
  Subroutine b_fun(a,b)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)    :: a(angular_momentum + 1 : angular_momentum + 1 + series_size)
  REAL*8, DIMENSION(:)    :: b( -angular_momentum - 1 : -angular_momentum + series_size)
  REAL*8                  :: pl
  REAL*8                  :: argument
  INTEGER                 :: i 
  INTEGER                 :: l_upper 
  pl=pl_fun()
  b( -angular_momentum - int_one ) = zero
  b( -angular_momentum ) = one
  argument = two * eta_in
  l_upper = -angular_momentum + series_size
  DO i = -angular_momentum + int_one , angular_momentum
     b(i) = argument * b(i-1)- b(i-2)
     b(i) = b(i) / ( (i - angular_momentum - int_one ) * ( i + angular_momentum ) )
  END DO
  b(angular_momentum + 1) = zero
  DO i = angular_momentum + int_two , l_upper
     b(i) = argument * b(i-1)- b(i-2) - ( i + i - int_one) * pl * a(i)
     b(i) = b(i) / ( ( i - angular_momentum - int_one) * ( i + angular_momentum ) )
  END DO
  END Subroutine b_fun
!====================================================================
!====================================================================
!deck series_coefficients_regular_whittaker_function
!***begin prologue     series_coefficients_regular_whittaker_function
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***description        series expansion coefficients for regular Whittaker function, M(a,b,z)
!                      
!***references          nbs mathematical handbook (abromowitz and stegun)
!***routines called
!***end prologue       series_coefficients_regular_whittaker_function
  Subroutine series_coefficients_regular_whittaker_function(a_w,b_w)
  IMPLICIT NONE
  REAL*8                  :: a_w
  REAL*8                  :: b_w
  REAL*8, DIMENSION(:)    :: whit_coef(0 : series_size)
  INTEGER                 :: i 
  whit_coef(0) = one
  whit_coef(1) = a_w / b_w
  DO i = 2 , series_size
     whit_coef(i) = whit_coef(i-1) * ( a_w + i - int_one ) / ( b_w + i - int_one)
  END DO
  END Subroutine series_coefficients_regular_whittaker_function
!====================================================================
!====================================================================
END MODULE Special_Functions_Module
