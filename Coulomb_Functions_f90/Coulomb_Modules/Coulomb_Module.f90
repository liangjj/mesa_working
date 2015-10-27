!*********************************************************************
                      MODULE Coulomb_Module
                      USE input_output
!
                      IMPLICIT NONE
!
!                     CONSTANTS
!
  COMPLEX*16                                         :: eye               = (0.d0,1.d0)
  REAL*8                                             :: pi                = 3.1415926535897932384626433832795028841968D0
  REAL*8                                             :: sqrt2             = 1.4142135623730950488016887242096980785696d0
  REAL*8                                             :: invsqrt2          = .7071067811865475244008443621048490392848d0
  REAL*8                                             :: eulerc            = .577215664901532860606512090082402431042D0
  REAL*8                                             :: tm30              = 1.0D-30
  REAL*8                                             :: tm16              = 1.0D-16
  REAL*8                                             :: sqrt_2_div_pi     = .79788456080286535587989211986876373695173d0
  REAL*8                                             :: zero              =  0.d0
  REAL*8                                             :: quarter           = .25d0
  REAL*8                                             :: half              = .50d0
  REAL*8                                             :: one               = 1.0d0
  REAL*8                                             :: two               = 2.0d0
  REAL*8                                             :: three             = 3.0d0
  REAL*8                                             :: four              = 4.0d0
  REAL*8                                             :: five              = 5.0d0
  REAL*8                                             :: six               = 6.0d0
  REAL*8                                             :: seven             = 7.0d0
  REAL*8                                             :: eight             = 8.0d0
  REAL*8                                             :: nine              = 9.0d0
  REAL*8                                             :: ten               = 10.d0
  REAL*8                                             :: ten_2             = 1.0D2
  REAL*8                                             :: ten_4             = 1.0D4
  REAL*8                                             :: ten_6             = 1.0D6
  REAL                                               :: zero_r4           = 0.0E0
  REAL                                               :: half_r4           = 0.5E0
  REAL                                               :: one_r4            = 1.0E0
  REAL                                               :: six_r4            = 6.0E0
  REAL                                               :: ten_r4            = 1.0E1
  REAL                                               :: rl35              = 3.5E1
  REAL                                               :: aloge             = 0.4342945E0
  REAL*8                                             :: abort             = 2.0D+04
  REAL*8                                             :: abort2            = 4.0D4
  INTEGER                                            :: int_one           = 1
  INTEGER                                            :: int_two           = 2
  INTEGER                                            :: int_three         = 3
  INTEGER                                            :: int_four          = 4
  INTEGER                                            :: int_five          = 5
  INTEGER                                            :: int_six           = 6
  INTEGER                                            :: int_seven         = 7
  INTEGER                                            :: int_eight         = 8
  INTEGER                                            :: int_nine          = 9
  INTEGER                                            :: int_ten           = 10
  INTEGER                                            :: iuo               = 70
  INTEGER                                            :: max_fact          = 100
!

!                     MAJOR VARIABLES
!
  REAL*8                                             :: energy
  REAL*8                                             :: charge
  REAL*8                                             :: k
  REAL*8                                             :: r_series 
  REAL*8                                             :: r_asymptotic
  LOGICAL                                            :: print_sigma_l
  LOGICAL                                            :: print_long_range_coefficients
  LOGICAL                                            :: print_convergence
  LOGICAL                                            :: print_short_range_coefficients
  INTEGER                                            :: series_size
  INTEGER                                            :: asymptotic_size
  REAL*8                                             :: eta_in
  REAL*8                                             :: wronskian
  CHARACTER(LEN=80)                                  :: title
  CHARACTER(LEN=1600)                                :: card
  CHARACTER(LEN=80)                                  :: cpass
  CHARACTER(LEN=80)                                  :: quantities_returned
  CHARACTER(LEN=80)                                  :: type
  INTEGER                                            :: l_val
  INTEGER                                            :: l_max
  INTEGER                                            :: l_min
  INTEGER                                            :: number_of_r_values
  INTEGER                                            :: iexp
  INTEGER                                            :: ifail
  REAL*8                                             :: r_min
  REAL*8                                             :: r_max
  REAL*8                                             :: r_step
  REAL*8                                             :: fl
  REAL*8                                             :: dfl
  REAL*8                                             :: gl
  REAL*8                                             :: dgl
  REAL*8                                             :: u
  REAL*8                                             :: twou
  REAL*8                                             :: fouru

!
!                    ALLOCATABLES
!
  REAL*8, DIMENSION(:,:), ALLOCATABLE                :: fgd
  REAL*8, DIMENSION(:),   ALLOCATABLE                :: rho
  REAL*8, DIMENSION(:),   ALLOCATABLE                :: rho_inv
  REAL*8, DIMENSION(:),   ALLOCATABLE                :: fact
  TYPE series_coefficients
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: a
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: b
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: a_0
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: b_0
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: c_0
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: d_0
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: a_i
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: b_i
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: e_0
  END TYPE series_coefficients
  TYPE(series_coefficients),                                         &
                  DIMENSION(:),                                      &
                  ALLOCATABLE                        :: power_series
  TYPE(series_coefficients),                                         &
                  DIMENSION(:),                                      &
                  ALLOCATABLE                        :: asymptotic_series
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  Subroutine factorials(fact)
  IMPLICIT NONE
  REAL*8, DIMENSION(0:max_fact)           :: fact
  INTEGER                                 :: i
  fact(0) = one
  DO i = 1, max_fact
     fact(i) = i * fact(i-1)
  END DO
  END SUBROUTINE factorials
!=======================================================================
!=======================================================================
!deck c0_sq_fn
  REAL*8 FUNCTION c0_sq_fn()
!***begin prologue     c0_sq_fn
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

!***end prologue       c0_sq_fn
  IMPLICIT NONE
  REAL*8 two_pi_eta
  two_pi_eta = two * pi * eta_in
  c0_sq_fn = two_pi_eta/(EXP(two_pi_eta)-one)
  END FUNCTION c0_sq_fn
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
  COMPLEX *16                          :: arg
  COMPLEX *16                          :: cgamma
  INTEGER                              :: l_1
  INTEGER                              :: two_el
  cl_fun = (two**l_val)*EXP(-pi*eta_in*half)
  l_1 = l_val + 1
  two_el = l_val + l_val
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
  dl_fun = one/( (l_val + l_val + 1)*cl_fun() )
  END FUNCTION dl_fun
!=======================================================================
!=======================================================================
!deck plfun
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
  REAL*8                   :: eta_sq
  INTEGER                  :: i
  INTEGER                  :: two_el
  INTEGER                  :: two_el_1
  INTEGER                  :: two_el_m_1
  pl_fun=two*eta_in
  IF (l_val > 0) THEN
      eta_sq = eta_in*eta_in
      pl_fun = pl_fun*(one+eta_sq)/three
      IF (l_val > 1) THEN
          DO i = 2, l_val
             two_el = i + i
             two_el_1 = two_el + 1
             two_el_m_1 = two_el-1
             pl_fun = pl_fun*four*(l_val*l_val + eta_sq)       &
                                /                              &
                      (two_el_1*two_el*two_el*two_el_m_1)
          END DO
      END IF
  END IF
  END FUNCTION pl_fun
!=======================================================================
!=======================================================================
!deck ql_pl_fn
!***begin prologue     ql_pl_fn
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***description        q /p  function for construction of irregular positive
!                       l  l energy coulomb functions
!***references          nbs mathematical handbook (abromowitz and stegun)
!***routines called    cpsi(complex psi function:clams)
!***end prologue       ql_pl_fn
  REAL*8 FUNCTION ql_pl_fn()
  IMPLICIT NONE
  COMPLEX*16                             :: cpsi
  COMPLEX*16                             :: eta_fac
  REAL*8                                 :: eta_sq
  REAL*8                                 :: sum_1
  REAL*8                                 :: sum_2
  INTEGER                                :: i 
  eta_fac = one + eta_in*eye
  IF (l_val == 0) THEN
      ql_pl_fn = - one + REAL( cpsi(eta_fac) ) + two*eulerc
  ELSE
      eta_sq = eta_in*eta_in
      sum_1 = zero
      DO i = 1 , l_val
         sum_1 = sum_1 + i/( i * i + eta_sq)
      END DO
      sum_2 = zero
      DO i = 1 ,l_val + l_val + 1
         sum_2 = sum_2 + one / i
      END DO
      ql_pl_fn = sum_1-sum_2 + REAL( cpsi(eta_fac) ) +two*eulerc +  &
                               rl_fun()/pl_fun()
  END IF
  END FUNCTION ql_pl_fn
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
  REAL*8                    :: denominator
  COMPLEX*16                :: com_eta
  COMPLEX*16                :: eta_fac
  COMPLEX*16                :: num
  COMPLEX*16                :: cr_l
  INTEGER                   :: i
  INTEGER                   :: n_terms
  INTEGER                   :: l_fac
  INTEGER                   :: i_test
  IF (l_val == 0) THEN
      rl_fun = zero
  ELSE
      com_eta = eta_in*eye
      n_terms = l_val+l_val+1
      l_fac = n_terms
      eta_fac = com_eta-l_val
      num = one
      cr_l = num/l_fac
      DO i= 2 , n_terms
         l_fac = l_fac-1
         num = num * two * eta_fac
         denominator = l_fac*fact(i-1)
         cr_l =cr_l + num / denominator
         eta_fac = eta_fac + one
      END DO
      rl_fun = imag(cr_l)
      i_test = l_val - 2 * (l_val/2)
      IF (i_test == 0) THEN
          rl_fun = -rl_fun
      END IF
      rl_fun = rl_fun/fact(l_val + l_val)
  END IF
  END FUNCTION rl_fun
!=======================================================================
!=======================================================================
!deck abcoef
  SUBROUTINE short_range_coefficients(a,b,a_0,b_0,c_0,d_0)
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
  REAL*8, DIMENSION(:)                     :: a(l_val+1 : l_val+1+series_size)
  REAL*8, DIMENSION(:)                     :: b(-l_val-1 : -l_val+series_size)
  REAL*8, DIMENSION(:)                     :: a_0(0 : series_size)
  REAL*8, DIMENSION(:)                     :: b_0(0 : series_size)
  REAL*8, DIMENSION(:)                     :: c_0(0 : series_size)
  REAL*8, DIMENSION(:)                     :: d_0(0 : series_size)
  REAL*8                                   :: gamma
  REAL*8                                   :: psi
  INTEGER                                  :: i
  INTEGER                                  :: l_1
  INTEGER                                  :: l_2
  INTEGER                                  :: l_3
  INTEGER                                  :: l_up_a
  INTEGER                                  :: l_up_b
  INTEGER                                  :: l_down
  INTEGER                                  :: ia
  INTEGER                                  :: ib
  INTEGER                                  :: two_el
  INTEGER                                  :: two_el_1
  INTEGER                                  :: two_el_2
  REAL*8                                   :: two_eta
  REAL*8                                   :: arg_dum
  REAL*8                                   :: pl
  CHARACTER(LEN=3)                         :: itoc
  l_1 = l_val + 1
  l_2 = l_val + 2
  l_3 = l_val + 3
  l_up_a = l_1 + series_size
  l_up_b = -l_val + series_size
  l_down = -l_val - 1
  IF (energy == 'positive') THEN
      two_eta = two*eta_in
!**********************************************************************c
!             generate the a coefficients needed for the regular       c
!                            series solution                           c
!**********************************************************************c
     a(l_1) = one
     a(l_2)= eta_in / l_1
     DO i = l_3 , l_up_a
        a(i) = two_eta * a(i-1) - a(i-2)
        a(i) = a(i)/ ( (i + l_val ) * ( i - l_1) )
     END DO
!**********************************************************************c
!              now for the b coefficients needed for the irregular     c
!              solution, things are more complicated since             c
!              they depend on the previous generation of the a's       c
!**********************************************************************c
     pl=pl_fun()
     b( -l_down ) = zero
     b(-l_val) = one
     DO i = -l_val+1 , l_val
        b(i) = two_eta * b(i-1)- b(i-2)
        b(i) = b(i) / ( (i - l_1 ) * ( i + l_val ) )
     END DO
     DO i = l_1 , -l_val+series_size
        b(i) = two_eta * b(i-1)- b(i-2) - ( i + i - 1) * pl * a(i)
        b(i) = b(i) / ( (i - l_1 ) * ( i + l_val ) )
     END DO
     IF (print_short_range_coefficients) THEN
         WRITE(iout,*)
         WRITE(iout,*) '          the a and b short range expansion '//      &
                       'coefficients for L = '//itoc(l_val)
         WRITE(iout,*)
         WRITE(iout,*) '      ia          a           ib'// '          b'
         DO i = 0, series_size
            ia = l_1 + i
            ib = -l_val + i
            WRITE(iout,1) ia, a(ia), ib,b(ib)
         END DO
     END IF
  ELSE IF (energy == 'negative') THEN
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
     arg_dum = eta_in + l_1
     a_0(0) = gamma(arg_dum)
     two_el = l_val +l_val 
     two_el_1 = two_el + 1
     two_el_2 = two_el + 2
     arg_dum = two_el_2
     a_0(0)=a_0(0)/gamma(arg_dum)
     DO i = 1 , series_size
        a_0(i) = ( i + eta_in + l_val) * a_0(i-1)/ ( i * ( i + two_el_1 ) )
     END DO
     arg_dum = 1.d0
     d_0(0) = psi(arg_dum)
     arg_dum = two_el_2
     d_0(0) = d_0(0) + psi(arg_dum)
     arg_dum = eta_in + l_1
     d_0(0)=d_0(0)- psi(arg_dum)
     DO i  = 1, series_size
        d_0(i) = d_0(i-1) + 1.d0/i + 1.d0/(i + two_el + 1.d0)               &
                                 - 1.d0/( i + eta_in + l_val)
     END DO
     DO i  = 1, series_size
        b_0(i) = a_0(i)*d_0(i)
     END DO
     arg_dum = two_el_1
     c_0(0)=gamma(eta_in-l_val)*gamma(arg_dum)
     IF (l_val > 0) THEN
         DO i = 1, l_1
            c_0(i)=( i + eta_in - l_1) * c_0(i-1) / ( i * (two_el_1-i) )
       END DO
     END IF
     IF (print_long_range_coefficients) THEN
         WRITE(iout,*)
         WRITE(iout,*) '          the a_0 , b_0 ,c_0 and d_0 long range '//  &
                       'expansion '// 'coefficients for L = '//itoc(l_val)
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
  END IF
1 FORMAT(5X,i3,3X,e15.8,3X,i3,3X,e15.8)
2 FORMAT(5X,i3,3X,e15.8,3X,e15.8)
3 FORMAT(5X,i3,21X,e15.8)
RETURN
END SUBROUTINE short_range_coefficients









!                                                                      *
!  C O U L F G  -  P A C K A G E   (FROM  A.R. B A R N E T T)          *
!                                                                      *
!  THIS PACKAGE IS USED TO CALCULATE REGULAR AND IRREGULAR             *
!  COULOMB - AND BESSEL - FUNCTIONS                                    *
!                                                                      *
SUBROUTINE coulfg(xx,eta_cfg,fc,gc,dfc,dgc)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                        
!  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD             
!                                                                        
!  A. R. BARNETT           MANCHESTER  MARCH   1981                      
!                                                                        
!  ORIGINAL PROGRAM 'RCWFN'      IN    CPC  8 (1974) 377-395             
!                 + 'RCWFF'      IN    CPC 11 (1976) 141-142             
!  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314             
!  THIS VERSION WRITTEN UP       IN    CPC 27 (1982) 147-166           
!                                                                      
!                                                                      
!                                                                      
!                                                                      
!  THIS CODE IS A MODIFIED VERSION OF THAT PUBLISHED BY BARNETT IN     
!  CPC 27 (1982) 147-166. IT HAS BEEN RE-CODED IN FORTRAN 77 AND       
!  THE ACCURACY IS DETERMINED BY A MACHINE DEPENDENT PARAMETER         
!  ( SEE BELOW UNDER 'ACCURACY' ).                                     
!                                                                      
!                                                                      
!  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA_cfg (INCLUDING 0), 
!   AND INTEGER LAMBDA.GE.0 FOR A RANGE OF LAMBDA VALUES:                
!   L_MIN TO L_MAX.                                                      
!   STARTING ARRAY ELEMENT IS M1 = L_MIN+1                              
!                                                                      C
!  IF 'MODE' = 1  GET F,G,F',G'                                        C
!            = 2  GET F,G                                              C
!            = 3  GET F                                                C
!  IF 'KFN'  = 0  REAL       COULOMB FUNCTIONS ARE RETURNED            C
!            = 1  SPHERICAL   BESSEL      "      "     "               C
!            = 2  CYLINDRICAL BESSEL      "      "     "               C
!                                                                      C
!  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT                          C
!                                                                      C
!  ACCURACY                                                            C
!  ========                                                            C
!                                                                      C
!                                                                      C
!  THE ACCURACY IS DETERMINED BY THE VARIABLE: ACCUR                   C
!  ACCUR IS SET TO: MAX(U,1.0D-16), WHERE U IS A MACHINE DEPENDENT     C
!  QUANTITY DETERMINED BY THE SUBROUTINE MACHIN. U IS A MEASURE OF     C
!  THE MACHINE ACCURACY.                                               C
!  THE USER MUST CALL MACHIN BEFORE THE FIRST CALL TO COULFG.          C
!  IF THE USER'S MACHINE ALLOWS MORE PRECISION THAN 1.0D-16 AND IF     C
!  A PRECISION BETTER THAN 1.0D-16 IS REQUIRED, THEN ALTER THE VALUE   C
!  OF 'TM16' IN THE PARAMETER STATEMENT BELOW.                         C
!  IN THE OSCILLATING REGION X.GE.XTURN, WHERE                         C
!  XTURN = ETA1+SQRT( ETA1**2+ L_MIN*( L_MIN+1) ), SOLUTIONS ARE       C
!  OBTAINED TO AN ACCURACY ACCUR. HOWEVER IF X IS SUFFICIENTLY         C
!  SMALLER THAN XTURN, SO THAT G.GT.1.0D6, THEN SOLUTIONS ARE          C
!  OBTAINED USING A JWKB APPROXIMATION AND THE RESULTS WILL BE MUCH    C
!  LESS ACCURATE ( IN GENERAL THE JWKB APPROXIMATION PROVIDES RESULTS  C
!  TO BETTER THAN 1% ). IF THE JWKB APPROXIMATION IS USED, A WARNING   C
!  MESSAGE IS PRINTED OUT FOR THE USER'S INFORMATION.                  C
!                                                                      C
!                                                                      C
!   OVERFLOW/UNDERFLOW                                                 C
!   ==================                                                 C
!                                                                      C
!   TO AVOID UNDERFLOW/OVERFLOW WHEN THE JWKB APPROXIMATION IS USED    C
!   THE USER MUST SET THE PARAMETER 'IUO' IN SUBROUTINE JWKB.          C
!                                                                      C
!                                                                      C
!   IEXP ON OUTPUT                                                     C
!   ==============                                                     C
!                                                                      C
!   IF IEXP .GT.1 ON OUTPUT, THEN SCALED RESULTS EXIST IN THE ARRAYS   C
!   FC,GC,DFC AND DGC. THE TRUE SOLUTIONS ARE FC*10**(-IEXP);          C
!   DFC*10**(-IEXP); GC*10**(IEXP); AND DGC*10**(IEXP).                C
!                                                                      C
!                                                                      C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
IMPLICIT NONE
REAL*8                                   :: xx
REAL*8, DIMENSION(1:l_max+1)             :: fc
REAL*8, DIMENSION(1:l_max+1)             :: gc
REAL*8, DIMENSION(1:l_max+1)             :: dfc
REAL*8, DIMENSION(1:l_max+1)             :: dgc
REAL*8                                   :: eta_cfg
INTEGER                                  :: mode1
INTEGER                                  :: kfn
REAL*8                                   :: x  
REAL*8                                   :: xl
REAL*8                                   :: xll
REAL*8                                   :: xll1
REAL*8                                   :: xlm
REAL*8                                   :: xi
REAL*8                                   :: eta
REAL*8                                   :: accur
REAL*8                                   :: acc
REAL*8                                   :: acc4
REAL*8                                   :: acch
REAL*8                                   :: lxtra 
REAL*8                                   :: fjwkb
REAL*8                                   :: gjwkb
REAL*8                                   :: paccq
REAL*8                                   :: gh
REAL*8                                   :: gh2
REAL*8                                   :: hll
REAL*8                                   :: sl
REAL*8                                   :: rl2
REAL*8                                   :: phi
REAL*8                                   :: phi10
REAL*8                                   :: e2mm1
REAL*8                                   :: f
REAL*8                                   :: fcl
REAL*8                                   :: gcl
REAL*8                                   :: gcl1
REAL*8                                   :: gpl
REAL*8                                   :: fcl1
REAL*8                                   :: fpl
REAL*8                                   :: fcm
REAL*8                                   :: w
REAL*8                                   :: gam
REAL*8                                   :: p
REAL*8                                   :: q
REAL*8                                   :: wi
REAL*8                                   :: ar
REAL*8                                   :: ai
REAL*8                                   :: br
REAL*8                                   :: bi
REAL*8                                   :: a
REAL*8                                   :: b
REAL*8                                   :: c
REAL*8                                   :: dr
REAL*8                                   :: di
REAL*8                                   :: dp
REAL*8                                   :: dq
REAL*8                                   :: rl
REAL*8                                   :: el
REAL*8                                   :: pk
REAL*8                                   :: px
REAL*8                                   :: ek
REAL*8                                   :: pk1
REAL*8                                   :: d
REAL*8                                   :: df
REAL*8                                   :: tk
REAL*8                                   :: alpha
REAL*8                                   :: beta
INTEGER                                  :: mode
INTEGER                                  :: npq
INTEGER                                  :: nfp
INTEGER                                  :: l
INTEGER                                  :: lp
INTEGER                                  :: m1
INTEGER                                  :: l1
INTEGER                                  :: maxl
LOGICAL                                  :: etane0
LOGICAL                                  :: xlturn
SAVE
!     INP AND IOUT SPECIFY THE INPUT AND OUTPUT UNITS
COMMON  / steed /  paccq,nfp,npq,m1
!     COMMON BLOCK IS FOR INFORMATION ONLY.  NOT REQUIRED IN CODE
!     COULFG HAS CALLS TO: ABS,ANINT,DBLE,MAX,MIN,SIGN,NINT,SQRT
!     CHECK DIMENSIONS
IF ( size(fc) < (l_max+1) ) THEN
  WRITE(iout,1007)
  WRITE(iout,1006) l_max
  STOP
END IF
ifail = 0
!     SET VALUES OF ACCURACY VARIABLES
accur = MAX(u,tm16)
acc = accur
acc4 = acc*ten_4
acch = SQRT(acc)
!    TEST RANGE OF XX, EXIT IF.LE.DSQRT(ACCUR) OR IF NEGATIVE
IF ( xx <= acch ) THEN
  ifail = -1
  WRITE(iout,*) 'This is IFAIL = -1 :'
  WRITE(iout,1007)
  WRITE(iout,1001) xx,acch
  RETURN
END IF
!     CHECK THAT L_MIN AND L_MAX ARE SUCH THAT:
!     XLM.GT.-ONE   AND   L_MAX.GE.L_MIN
xlm = ANINT(DBLE(l_min))
IF ( xlm <= -one .OR. l_max < l_min ) THEN
  ifail = -2
  WRITE(iout,*) 'This is IFAIL = -2 :'
  WRITE(iout,1007)
  WRITE(iout,1002) l_max,l_min,xlm
  RETURN
END IF
IF ( type == 'coulomb') THEN
     kfn = 0
ELSE IF (type == 'spherical_bessel') THEN
     kfn = 1
ELSE IF (type == 'cylindrical bessel') THEN
     kfn = 2
ELSE
     write(iout,100)
     stop
END IF
IF ( quantities_returned == 'functions_and_derivatives') THEN
     mode1 = 1
ELSE IF (quantities_returned == 'functions') THEN
     mode1 = 2
ELSE IF (quantities_returned == 'regular_function') THEN
     mode1 = 3
ELSE
     write(iout,200)
     stop
END IF
IF ( kfn == 2 ) THEN
  xlm = xlm-half
END IF
!     DETERMINE LXTRA = THE NUMBER OF ADDITIONAL LAMBDA VALUES
!     TO BE COMPUTED.
lxtra = l_max-l_min

IF ( mode1 == 1 ) THEN
  mode = 1
ELSE
  mode = mode1
END IF
iexp = 1
npq = 0
gjwkb = zero
paccq = one
IF ( kfn /= 0 ) THEN
  eta = zero
  etane0 = .false.
ELSE
  eta = eta_cfg
  etane0 = .true.
END IF
x = xx
!     DETERMINE WHETHER X IS .LT. THE TURNING POINT VALUE
!     IF IT IS: SET XLTURN = .TRUE.*     IF NOT  : SET XLTURN = .FALSE.
IF ( x*(x - two*eta) < xlm*xlm + xlm ) THEN
  xlturn = .true.
ELSE
  xlturn = .false.
END IF
!!!!      write(*,*) 'XLTURN =', XLTURN
e2mm1 = eta*eta + xlm*xlm + xlm
xll = xlm + DBLE(lxtra)
!       XLL  ISMAX LAMBDA VALUE, OR 0.5 SMALLER FOR J,Y BESSELS
!         DETERMINE STARTING ARRAY ELEMENT (M1) FROM L_MIN
m1 = l_min+1
l1 = m1 + lxtra
!    EVALUATE CF1  = F    =  FPRIME(XL,ETA,X)/F(XL,ETA,X)
xi = one/x
fcl = one
pk = xll + one
px = pk + abort
1     CONTINUE
ek = eta / pk
f = (ek + pk*xi)*fcl + (fcl - one)*xi
pk1 = pk + one
!   ENSURE THAT B1 .NE. ZERO FOR NEGATIVE ETA: FIXUP IS EXACT.
IF ( ABS(eta*x + pk*pk1) > acc ) THEN
  GO TO 2
END IF
fcl = (one + ek*ek)/(one + (eta/pk1)**2)
pk = two + pk
GO TO 1
2     d = one/((pk + pk1)*(xi + ek/pk1))
df = -fcl*(one + ek*ek)*d
IF ( fcl /= one ) THEN
  fcl = -one
END IF
IF ( d < zero) THEN
  fcl = -fcl
END IF
f = f + df
!   BEGIN CF1 LOOP ON PK = K = LAMBDA + 1
p = one
3     CONTINUE
pk = pk1
pk1 = pk1 + one
ek = eta / pk
tk = (pk + pk1)*(xi + ek/pk1)
d = tk - d*(one + ek*ek)
!!!          WRITE(IOUT,1003) ABORT,F,DF,PK,PX,ACC
!!!          STOP
IF ( ABS(d) <= acch ) THEN
  WRITE (iout,1000) d,df,acch,pk,ek,eta,x
  p = p + one
  IF( p > two ) THEN
    ifail = 1
    WRITE(iout,*) 'This is IFAIL = 1 :'
    WRITE(iout,1007)
    WRITE(iout,1003) abort,f,df,pk,px,acc
    RETURN
  END IF
END IF
d = one/d
IF ( d < zero ) THEN
  fcl = -fcl
END IF
df = df*(d*tk - one)
f = f + df
IF ( pk > px ) THEN
  ifail = 1
  WRITE(iout,*) 'This is IFAIL = 1 :'
  WRITE(iout,1007)
  WRITE(iout,1003) abort,f,df,pk,px,acc
  RETURN
END IF
IF ( ABS(df) < ABS(f)*acc ) THEN
  GO TO 4
END IF
GO TO 3
4     CONTINUE
nfp = nint(pk - xll - one)
IF ( lxtra /= 0 ) THEN
!   DOWNWARD RECURRENCE TO LAMBDA = XLM. ARRAY GC,IF PRESENT,STORES RL
  fcl = fcl*tm30
  fpl = fcl*f
  IF ( mode == 1 ) THEN
    dfc(l1) = fpl
  END IF
  fc (l1) = fcl
  xl = xll
  rl = one
  el = zero
  DO   lp = 1,lxtra
    IF ( etane0 ) THEN
      el = eta/xl
      rl = SQRT(one+el*el)
    END IF
    sl = el + xl*xi
    l = l1 - lp
    fcl1 = (fcl *sl + fpl)/rl
    fpl = fcl1*sl - fcl *rl
    fcl = fcl1
    fc(l) = fcl
    IF ( mode == 1 ) THEN
      dfc(l) = fpl
    END IF
    IF ( mode /= 3 .AND. etane0 ) THEN
      gc(l+1) = rl
    END IF
    xl = xl - one
  END DO
  IF ( fcl == zero ) THEN
    fcl = acc
  END IF
  f = fpl/fcl
END IF
!    NOW WE HAVE REACHED LAMBDA = XL_MIN = XLM
!    EVALUATE CF2 = P + I.Q  AGAIN USING STEED'S ALGORITHM
!    SEE TEXT FOR COMPACT COMPLEX CODE FOR SP CDC OR NON-ANSI IBM
IF ( xlturn ) THEN
  CALL jwkb(x,eta,MAX(xlm,zero),fjwkb,gjwkb,iexp)
END IF
IF( iexp > 1 .OR. gjwkb > ten_6 ) THEN
!     AT THIS POINT  G(XLM) .GT. 10**6 OR IEXP .GT. IUO & XLTURN = .TRUE
  WRITE(iout,1008)
  w = fjwkb
  gam = gjwkb*w
  p = f
  q = one
ELSE
  xlturn = .false.
  pk = zero
  wi = eta + eta
  p = zero
  q = one - eta*xi
  ar = -e2mm1
  ai = eta
  br = two*(x - eta)
  bi = two
  dr = br/(br*br + bi*bi)
  di = -bi/(br*br + bi*bi)
  dp = -xi*(ar*di + ai*dr)
  dq = xi*(ar*dr - ai*di)
  6       CONTINUE
  p = p + dp
  q = q + dq
  pk = pk + two
  ar = ar + pk
  ai = ai + wi
  bi = bi + two
  d = ar*dr - ai*di + br
  di = ai*dr + ar*di + bi
  c = one/(d*d + di*di)
  dr = c*d
  di = -c*di
  a = br*dr - bi*di - one
  b = bi*dr + br*di
  c = dp*a - dq*b
  dq = dp*b + dq*a
  dp = c
  IF ( pk > abort2 ) THEN
    ifail = 2
    WRITE(iout,*) 'This is IFAIL = 2 :'
    WRITE(iout,1007)
    WRITE(iout,1004) abort,p,q,dp,dq,acc
    RETURN
  END IF
  IF ( ABS(dp)+ABS(dq) < (ABS(p)+ABS(q))*acc ) THEN
    GO TO 7
  END IF
  GO TO 6
  7       CONTINUE
  npq = nint(pk/two)
  paccq = half*acc/MIN(ABS(q),one)
  IF ( ABS(p) > ABS(q) ) THEN
    paccq = paccq*ABS(p)
  END IF
  gam = (f - p)/q
  IF ( q <= acc4*ABS(p) ) THEN
    ifail = 3
    WRITE(iout,*) 'This is IFAIL = 3 :'
    WRITE(iout,1007)
    WRITE(iout,1005) p,q,acc,lxtra,m1
    RETURN
  END IF
  w = one/SQRT((f - p)*gam + q)
END IF
!    NORMALISE FOR SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
IF ( kfn == 0 ) THEN
  alpha = zero
  beta = one
ELSE  &
  IF ( kfn == 1 ) THEN
    alpha = xi
    beta = xi
  ELSE
    alpha = xi*half
    beta = SQRT(xi)*sqrt_2_div_pi
  END IF
  fcm = SIGN(w,fcl)*beta
  fc(m1) = fcm
  IF ( mode < 3 ) THEN
    IF ( xlturn ) THEN
      gcl = gjwkb*beta
    ELSE
      gcl = fcm*gam
    END IF
    IF ( kfn /= 0 ) THEN
      gcl = -gcl
    END IF
    gc(m1) = gcl
    gpl = gcl*(p - q/gam) - alpha*gcl
    IF ( mode == 1 ) THEN
      dgc(m1) = gpl
      dfc(m1) = fcm*(f - alpha)
    END IF
  END IF
  IF ( lxtra /= 0 ) THEN
!     UPWARD RECURRENCE FROM GC(M1),DGC(M1)  STORED VALUE IS RL
!     RENORMALISE FC,DFC AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
!        XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
    w = beta*w/ABS(fcl)
    maxl = l1 - 1
    DO   l = m1,maxl
      IF ( mode < 3 ) THEN
        xl = xl + one
        IF ( etane0 ) THEN
          el = eta/xl
          rl = gc(l+1)
        END IF
        sl = el + xl*xi
        gcl1 = ((sl - alpha)*gcl - gpl)/rl
        gpl = rl*gcl - (sl + alpha)*gcl1
        gcl = gcl1
        gc(l+1) = gcl1
        IF ( mode == 1 ) THEN
          dgc(l+1) = gpl
          dfc(l+1) = w*(dfc(l+1) - alpha*fc(l+1))
        END IF
      END IF
      fc(l+1) = w* fc(l+1)
    END DO
  END IF
  RETURN
!     FORMAT STATEMENTS:
  1000  FORMAT (/' CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ',7D 9.2/)
  1001  FORMAT(' FOR XX = ',d12.3,' TRY SMALL-X  SOLUTIONS',  &
      ' OR X NEGATIVE'/ ,' SQUARE ROOT ACCURACY PARAMETER = ',d12.3/)
  1002  FORMAT(/' PROBLEM WITH INPUT ORDER VALUES:L_MAX,L_MIN,XLM = ',  &
      2I6,d15.6/)
  1003  FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',f10.0,' ITERATIONS',/  &
      ' F,DF,PK,PX,ACCUR =  ',5D12.3//)
  1004  FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',f7.0,' ITERATIONS',/  &
      ' P,Q,DP,DQ,ACCUR =  ',4D17.7,d12.3//)
  1005  FORMAT(' FINAL Q.LE.DABS(P)*ACC*10**4 , P,Q,ACC = ',3D12.3,4X,  &
      ' LXTRA,M1 = ',2I5 /)
  1006  FORMAT(' ARRAY DIMENSIONS SHOULD BE INCREASED TO INT(L_MAX+1)')
  1007  FORMAT(//' COULFG FAILURE. '//)
  1008  FORMAT(//' INFORMATION MESSAGE.'/  &
      ' SOLUTIONS WERE OBTAINED USING THE JWKB APPROXIMATION '/  &
      ' THEY MAY BE LESS ACCURATE THAN DESIRED'//)
  100   FORMAT(//'FUNCTION DESIGNATION INCORRECT')
  200   FORMAT(//'INCORRECT REQUEST')
END SUBROUTINE coulfg
! ==============================================================
! ==============================================================
SUBROUTINE jwkb(xx,eta1,xl,fjwkb,gjwkb,iexp)
IMPLICIT NONE
REAL*8                         :: xx
REAL*8                         :: eta1
REAL*8                         :: xl
REAL*8                         :: fjwkb
REAL*8                         :: gjwkb
INTEGER                        :: iexp
REAL                           :: eta
REAL                           :: x
REAL                           :: gh2
REAL                           :: xll1
REAL                           :: hll
REAL                           :: hl
REAL                           :: sl
REAL                           :: rl2
REAL                           :: gh
REAL                           :: phi
REAL                           :: phi10


!     COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS   FOR XL.GE. 0
!     AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
!     CALLS ATAN2,EXP,INT,LOG,MAX,REAL,SQRT

!    IUO MUST BE SET BY THE USER. ITS PURPOSE IS TO PREVENT
!    UNDERFLOW/OVERFLOW AND SO IT IS MACHINE DEPENDENT.
!    IUO SHOULD HAVE AS VALUE THE EXPONENT OF THE LARGEST
!    REAL NUMBER WHICH THE MACHINE IS CAPABLE OF HOLDING
!    WITHOUT OVERFLOW MINUS 5.


x     = xx
eta   = eta1
gh2   = x*(eta + eta - x)
xll1  = MAX(xl*xl + xl,zero)
IF ( gh2 + xll1 > zero_r4 ) THEN
  hll = xll1 + six_r4/rl35
  hl = SQRT(hll)
  sl = eta/hl + hl/x
  rl2 = one_r4 + eta*eta/hll
  gh = SQRT(gh2 + hll)/x
  phi = x*gh - half_r4*( hl*LOG((gh + sl)**2/rl2) - LOG(gh) )
  IF ( eta /= zero_r4 ) THEN
    phi = phi - eta*ATAN2(x*gh,x - eta)
  END IF
  phi10 = -phi*aloge
  iexp = INT(phi10)
  IF ( iexp > iuo ) THEN
    gjwkb = ten_r4**(phi10 - REAL(iexp))
  ELSE
    gjwkb = EXP(-phi)
    iexp = 0
  END IF
  fjwkb = half_r4/(gh*gjwkb)
END IF
RETURN
END SUBROUTINE jwkb
!====================================================================
!====================================================================
SUBROUTINE machin
!     U IS THE SMALLEST POSITIVE NUMBER SUCH THAT
!     (1.+U) .GT. 1.
!     U IS COMPUTED APPROXIMATELY AS A POWER OF 1./2.
!     THIS ROUTINE IS COMPLETELY EXPLAINED AND DOCUMENTED
!     IN THE TEXT:
!     " COMPUTER SOLUTION OF ORDINARY DIFFERENTIAL EQUATIONS:
!       THE INITIAL VALUE PROBLEM " BY L.F.SHAMPINE AND M.K.GORDON
!     THE ROUTINE HAS BEEN RE-CODED IN FORTRAN 77
IMPLICIT DOUBLE PRECISION (a-h,o-z)
SAVE
halfu = half
1     t = one+halfu
IF ( t <= one ) THEN
  u = two*halfu
  twou = two*u
  fouru = two*twou
ELSE
  halfu = half*halfu
  GO TO 1
END IF
RETURN
END SUBROUTINE machin
!====================================================================
!====================================================================
END MODULE Coulomb_Module
