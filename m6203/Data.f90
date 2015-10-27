!deck Data
!***begin prologue     Data
!***date written       140601   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Data
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            Data for the numerical grid generation code
!***description  
!***             
!***             
!***references
!***routines called    
!***end prologue       Data
!********************************************************************************
!********************************************************************************
                        MODULE Data
  USE accuracy
  USE input_output,     ONLY: inp, iout
  USE Lebedev_Data
  IMPLICIT NONE
!********************************************************************************
!********************************************************************************
   REAL(idp)                              ::          pi    = 3.141592653589793238462643383276D0 
   REAL(idp)                              ::          two_pi= 6.283185307179586476925286766552D0 
   REAL(idp)                              ::          four_pi= 12.566370614359172953850573533104D0 
   REAL(idp)                              ::          zero    =  0.D0
   REAL(idp)                              ::          quarter = .25D0
   REAL(idp)                              ::          half    = .50D0
   REAL(idp)                              ::          third   = .33333333333333333333333333333D0
   REAL(idp)                              ::          fourth  = .25000000000000000000000000000D0
   REAL(idp)                              ::          fifth   = .20000000000000000000000000000D0
   REAL(idp)                              ::          sixth   = .16666666666666666666666666666D0
   REAL(idp)                              ::          seventh = .14285714285714285714285714285D0
   REAL(idp)                              ::          eighth  = .12500000000000000000000000000D0
   REAL(idp)                              ::          ninth   = .11111111111111111111111111111D0
   REAL(idp)                              ::          tenth   = .10000000000000000000000000000D0
   REAL(idp)                              ::          one     = 1.0D0
   REAL(idp)                              ::          two     = 2.0D0
   REAL(idp)                              ::          three   = 3.0D0
   REAL(idp)                              ::          four    = 4.0D0
   REAL(idp)                              ::          five    = 5.0D0
   REAL(idp)                              ::          six     = 6.0D0
   REAL(idp)                              ::          seven   = 7.0D0
   REAL(idp)                              ::          eight   = 8.0D0
   REAL(idp)                              ::          nine    = 9.0D0
   REAL(idp)                              ::          ten     = 10.D0
   INTEGER                                ::          int_zero   = 0
   INTEGER                                ::          int_one    = 1
   INTEGER                                ::          int_two    = 2
   INTEGER                                ::          int_three  = 3
   INTEGER                                ::          int_four   = 4
   INTEGER                                ::          int_five   = 5
   INTEGER                                ::          int_six    = 6
   INTEGER                                ::          int_seven  = 7
   INTEGER                                ::          int_eight  = 8
   INTEGER                                ::          int_nine   = 9
   INTEGER                                ::          int_ten    = 10
   INTEGER                                ::          int_eleven    = 11
   INTEGER                                ::          int_twelve    = 12
   INTEGER                                ::          int_thirteen  = 13
   INTEGER                                ::          int_fourteen  = 14
   INTEGER                                ::          int_fifteen   = 15
   INTEGER                                ::          int_sixteen   = 16
   INTEGER                                ::          int_seventeen = 17
   INTEGER                                ::          int_eighteen  = 18
   INTEGER                                ::          int_nineteen  = 19
   INTEGER                                ::          int_twenty    = 20
   REAL(idp)                              ::          nrzero  = 1.D-06
!
!                    hbar in joule-sec
!
   REAL(idp)                              ::          hbar = 1.054571596D-34
!
!                    electron mass in kg
!
   REAL(idp)                              ::          massau = 9.10938188D-31
!
!                    bohr radius in meters
!
   REAL(idp)                              ::          lenau = 5.291772083D-11
!
!                    time for an electron to make one bohr orbit in seconds
!
   REAL(idp)                              ::          timau    = 2.418884326D-17
   REAL(idp)                              ::          efieldau = 5.14220624D+11
   REAL(idp)                              ::          electric_field_to_intensity = 3.509338D+16
   REAL(idp)                              ::          peak_electric_field = .2849540283D-03
   REAL(idp)                              ::          pmass    = 1.67262158D-27
   REAL(idp)                              ::          massn2p  = 1.00137841887D0
   REAL(idp)                              ::          au_in_ev = 27.211396132D0
  CHARACTER (LEN=4096)                    :: ops
  CHARACTER (LEN=80)                      :: cpass
  CHARACTER (LEN=30)                      :: str
  CHARACTER (LEN=128)                     :: grid_filename
  CHARACTER (LEN=3)                       :: chra
  CHARACTER (LEN=3)                       :: yn
  CHARACTER (LEN=1600)                    :: card
  CHARACTER (LEN=8)                       :: rad_grid
  CHARACTER (LEN=8)                       :: ang_grid
  REAL(idp)                               :: alpha
  REAL(idp)                               :: cutoff
  REAL(idp)                               :: depth
  REAL(idp)                               :: charge
  REAL(idp)                               :: awell
  REAL(idp)                               :: e_c
  REAL(idp)                               :: mass
  REAL(idp)                               :: factor
  REAL(idp)                               :: omega
  REAL(idp)                               :: shift
  REAL(idp), DIMENSION(2)                 :: amp
  REAL(idp), DIMENSION(2)                 :: expnt
  REAL(idp)                               :: pre_factor
  REAL(idp)                               :: dscale
  REAL(idp)                               :: exact_vol
  REAL(idp)                               :: yukawa_integral
  INTEGER                                 :: n_p
  INTEGER                                 :: n_scale
  INTEGER                                 :: nwell
  INTEGER                                 :: nocen    = 300
  INTEGER                                 :: numshl   = 1000
  INTEGER                                 :: no_radial_reg
  INTEGER                                 :: no_theta_reg
  INTEGER                                 :: no_phi_reg
  INTEGER                                 :: nedges
  INTEGER                                 :: ncent
  INTEGER                                 :: ncplus
  INTEGER                                 :: niter
  INTEGER                                 :: nrmax
  INTEGER                                 :: nthmax
  INTEGER                                 :: nphmax
  INTEGER                                 :: maxgrd
  INTEGER                                 :: maxshl
  INTEGER                                 :: ntotal
  INTEGER                                 :: nwttot
  INTEGER                                 :: ltop
  INTEGER                                 :: mtop
  INTEGER                                 :: lower
  INTEGER                                 :: upper
  INTEGER                                 :: skip
  INTEGER                                 :: maxang
  INTEGER                                 :: lstrng
  LOGICAL                                 :: no_voroni
  LOGICAL                                 :: no_scat
  LOGICAL                                 :: nonsep
  LOGICAL                                 :: yukawa_on
  LOGICAL                                 :: no_disk
  LOGICAL                                 :: one_grid
  LOGICAL                                 :: make_3d_grid
  LOGICAL                                 :: fixed_angular_quadrature
  LOGICAL, DIMENSION(10)                  :: print_leb
  LOGICAL, DIMENSION(10)                  :: print_ang
  LOGICAL, DIMENSION(10)                  :: print_gauss
  LOGICAL, DIMENSION(10)                  :: print_norm
  LOGICAL, DIMENSION(10)                  :: print_rad
  LOGICAL, DIMENSION(10)                  :: print_mat
  LOGICAL                                 :: print_yukawa
  LOGICAL                                 :: fedvr_functions
  LOGICAL                                 :: iosys_on
!********************************************************************************
!********************************************************************************
  END MODULE DATA
!********************************************************************************
!********************************************************************************
