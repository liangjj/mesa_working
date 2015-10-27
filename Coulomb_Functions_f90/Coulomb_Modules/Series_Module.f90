!*********************************************************************
                      MODULE Series_Module
                      USE Coulomb_Variables_Module
                      USE Special_Functions_Module
!
                      IMPLICIT NONE
  INTEGER                                  :: l_2
  INTEGER                                  :: l_3
  INTEGER                                  :: l_up_a
  INTEGER                                  :: l_up_b
  INTEGER                                  :: l_down
  INTEGER                                  :: ia
  INTEGER                                  :: ib
  REAL*8                                   :: two_eta
  REAL*8                                   :: pl
  REAL*8, EXTERNAL                         :: psi
  INTEGER                                  :: two_el_2
  REAL*8                                   :: arg_dum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        Contains
!=======================================================================
!=======================================================================
!deck positive_energy_short_range_coefficients
  SUBROUTINE positive_energy_short_range_coefficients(a,b)
!***begin prologue     short_range_coefficients
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            expansion coefficients for coulomb functions
!***                   at small rho.
!***description
!***references         nbs mathematical handbook
!***routines called
!***end prologue       short_range_coefficients
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                     :: a(angular_momentum+1 : angular_momentum+1+series_size)
  REAL*8, DIMENSION(:)                     :: b(-angular_momentum-1 : -angular_momentum+series_size)
  INTEGER                                  :: i  
!**********************************************************************c
!             generate the a coefficients needed for the regular       c
!                            series solution                           c
!**********************************************************************c
   Call a_fun(a)
!**********************************************************************c
!              now for the b coefficients needed for the irregular     c
!              solution, things are more complicated since             c
!              they depend on the previous generation of the a's       c
!**********************************************************************c
  Call b_fun(a,b)
  IF (print_short_range_coefficients) THEN
      WRITE(iout,*)
      WRITE(iout,*) '          the a and b short range expansion '//      &
                    'coefficients for L = '//itoc(angular_momentum)
      WRITE(iout,*)
      WRITE(iout,*) '     ia          a           ib'// '          b'
      DO i = 0, series_size
         ia = angular_momentum + int_one + i
         ib = -angular_momentum + i
         WRITE(iout,1) ia, a(ia), ib,b(ib)
      END DO
  END IF
1 FORMAT(5X,i3,3X,e15.8,3X,i3,3X,e15.8)
2 FORMAT(5X,i3,3X,e15.8,3X,e15.8)
3 FORMAT(5X,i3,21X,e15.8)
END SUBROUTINE positive_energy_short_range_coefficients
!=======================================================================
!=======================================================================
!deck negative_energy_short_range_coefficients
  SUBROUTINE negative_energy_short_range_coefficients(a_0,b_0,c_0,d_0)
!***begin prologue     short_range_coefficients
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            expansion coefficients for coulomb functions
!***                   at small rho.
!***description
!***references         nbs mathematical handbook
!***routines called
!***end prologue       short_range_coefficients
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                     :: a_0(0 : series_size)
  REAL*8, DIMENSION(:)                     :: b_0(0 : series_size)
  REAL*8, DIMENSION(:)                     :: c_0(0 : series_size)
  REAL*8, DIMENSION(:)                     :: d_0(0 : series_size)
  INTEGER                                  :: i  
!
!**********************************************************************c
!             coefficients for series solution for negative energy     c
!             regular and irregular coulomb functions at small         c
!                               distances.                             c
!          * the indexing is less fancy here and follows the paper     c
!                     of Henry and Roundtree in CPC                    c
!          * there are errors in that paper which will ne noted        c
!          * the coefficients defined here are identical to the        c
!            paper cited and the errors have been fixed in the         c
!            definition of the small rho and large rho behavior        c
!                           of the functions                           c
!**********************************************************************c
  l_1 = angular_momentum + 1
  arg_dum = eta_in + l_1
  a_0(0) = gamma(arg_dum)
  two_el = angular_momentum +angular_momentum 
  two_el_1 = two_el + 1
  two_el_2 = two_el + 2
  arg_dum = two_el_2
  a_0(0)=a_0(0)/gamma(arg_dum)
  DO i = 1 , series_size
     a_0(i) = ( i + eta_in + angular_momentum) * a_0(i-1)/ ( i * ( i + two_el_1 ) )
  END DO
  arg_dum = 1.d0
  d_0(0) = psi(arg_dum)
  arg_dum = two_el_2
  d_0(0) = d_0(0) + psi(arg_dum)
  arg_dum = eta_in + l_1
  d_0(0)=d_0(0)- psi(arg_dum)
  DO i  = 1, series_size
     d_0(i) = d_0(i-1) + 1.d0/i + 1.d0/(i + two_el + 1.d0)               &
                       - 1.d0/( i + eta_in + angular_momentum)
  END DO
  DO i  = 1, series_size
     b_0(i) = a_0(i)*d_0(i)
  END DO
  arg_dum = two_el_1
  c_0(0)=gamma(eta_in-angular_momentum)*gamma(arg_dum)
  IF (angular_momentum > 0) THEN
      DO i = 1, l_1
         c_0(i)=( i + eta_in - l_1) * c_0(i-1) / ( i * (two_el_1-i) )
    END DO
  END IF
  IF (print_short_range_coefficients) THEN
      WRITE(iout,*)
      WRITE(iout,*) '          the a_0 , b_0 ,c_0 and d_0 short range '//  &
                    'expansion '// 'coefficients for L = '//itoc(angular_momentum)
      WRITE(iout,*)
      WRITE(iout,*) '      i          a_0              b_0'
      DO i = 0 , series_size
         WRITE(iout,2) i,a_0(i),b_0(i)
      END DO
      WRITE(iout,*)
      WRITE(iout,*) '      i          c_0              d_0'
      DO i = 0 , series_size
         IF (i <= two_el) THEN
             WRITE(iout,2) i, c_0(i), d_0(i)
         ELSE
             WRITE(iout,3) i, d_0(i)
         END IF
      END DO
  END IF
1 FORMAT(5X,i3,3X,e15.8,3X,i3,3X,e15.8)
2 FORMAT(5X,i3,3X,e15.8,3X,e15.8)
3 FORMAT(5X,i3,21X,e15.8)
RETURN
END SUBROUTINE negative_energy_short_range_coefficients
END MODULE Series_Module
