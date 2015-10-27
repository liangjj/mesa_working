!*********************************************************************
                      MODULE Coulomb_Functions_Module
                      USE Coulomb_Variables_Module
                      USE Series_Module
!
                      IMPLICIT NONE
  REAL*8                                    :: cl
  REAL*8                                    :: pre
  REAL*8                                    :: add
  REAL*8                                    :: two_rho
  REAL*8                                    :: pre_fac
  REAL*8                                    :: dl
  REAL*8                                    :: argm
  REAL*8                                    :: exp_fac
  REAL*8                                    :: dxp_fac
  REAL*8                                    :: fac
  REAL*8                                    :: i_fac
  REAL*8                                    :: lg_fac
  REAL*8                                    :: fac_1
  REAL*8                                    :: fac_2
  REAL*8                                    :: fac_3
  REAL*8                                    :: c_sum
  REAL*8                                    :: dc_sum
  REAL*8                                    :: a_sum
  REAL*8                                    :: da_sum
  REAL*8                                    :: b_sum
  REAL*8                                    :: db_sum
  REAL*8                                    :: sigma_l
  REAL*8, DIMENSION(2)                      :: s_n
  REAL*8, DIMENSION(2)                      :: c_n
  REAL*8                                    :: s_sigma
  REAL*8                                    :: c_sigma
  REAL*8                                    :: theta_l
  REAL*8                                    :: theta
  REAL*8                                    :: f
  REAL*8                                    :: g
  REAL*8                                    :: f_d
  REAL*8                                    :: g_d
  REAL*8                                    :: f_old
  REAL*8                                    :: g_old
  REAL*8                                    :: f_new
  REAL*8                                    :: g_new
  REAL*8                                    :: f_old_d
  REAL*8                                    :: g_old_d
  REAL*8                                    :: f_d_new
  REAL*8                                    :: g_d_new
  INTEGER                                   :: count
  INTEGER                                   :: lower
  INTEGER                                   :: upper
  INTEGER                                   :: i_max
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        Contains
!=======================================================================
!=======================================================================
!deck series_expansion_regular_positive_energy_function
!***begin prologue     series_expansion_regular_positive_energy_function
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            series expansion for positive energy
!***                   regular coulomb function at small rho.
!***description
!***references         NBS handbook
!***routines called
!***end prologue       series_expansion_regular_positive_energy_function
  SUBROUTINE series_expansion_regular_positive_energy_function(a)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)        :: a(angular_momentum+1 : angular_momentum+1+series_size)
  INTEGER                     :: i
  pre = one
  count = 0
  lower = angular_momentum + int_one
  upper = angular_momentum + series_size
  fl = zero
  dfl = zero
  DO i = lower , upper
     count = count + int_one
     add = a(i)*pre
!**********************************************************************c
!            this funny business is done because accidental zeros      c
!            can occur in the a coefficients and at least a            c
!            reasonable result is likely to come out with 25 terms     c
!**********************************************************************c
     IF (ABS(add) < tol.AND.count >= 25) EXIT
         fl = fl + add
         dfl = dfl + i*add
         pre = pre*r
  END DO
  IF (print_convergence) THEN
      WRITE (iout,*) 'Power Series Converged For Regular Positive Energy Solution '//  &
                     'L = '//itoc(angular_momentum)//' in ', count,' Terms'
  END IF
  pre = r**angular_momentum
  cl = cl_fun()
  fl =  cl * fl * pre * r
  dfl = cl * dfl * pre
END SUBROUTINE series_expansion_regular_positive_energy_function
!=======================================================================
!=======================================================================
!deck series_expansion_irregular_positive_energy_function
!***begin prologue     series_expansion_irregular_positive_energy_function
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            series expansion for positive energy
!***                   irregular coulomb function at small rho.
!***description
!***references         NBS handbook
!***routines called
!***end prologue       series_expansion_irregular_positive_energy_function
  SUBROUTINE series_expansion_irregular_positive_energy_function(a,b,computed)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)        :: a(angular_momentum+1 : angular_momentum+1+series_size)
  REAL*8, DIMENSION(:)        :: b(-angular_momentum-1 : -angular_momentum+series_size)
  INTEGER                     :: i
  LOGICAL                     :: computed
!
!**********************************************************************c
!        gl = 2. * eta * fl * { ln(2.*rho) + ql/pl } / c0**2           c
!                      -dl * rho**(-l) * sum (n=0 to n=infinity)       c
!                                        b(n) * rho**n                 c
!**********************************************************************c
!**********************************************************************c
!        We need the regular solution to compute the irregular one     c
!
  IF (.not.computed) THEN
       Call series_expansion_regular_positive_energy_function(a)
  END IF
!
!**********************************************************************c
!**********************************************************************c
!                  do the theta series sum                             c
!**********************************************************************c
  pre = one
  count = 0
  gl = zero
  dgl = zero
  DO i = -angular_momentum, -angular_momentum + series_size
     count = count + int_one
     add = b(i)*pre
!**********************************************************************c
!            this funny business is done because accidental zeros      c
!            can occur in the a coefficients and at least a            c
!            reasonable result is likely to come out with 25 terms     c
!**********************************************************************c
    IF (ABS(add) < tol.AND.count >= 25) EXIT
    gl = gl + add
    dgl = dgl + i*add
    pre = pre*r
    IF (print_convergence) THEN
      WRITE (iout,*) 'Power Series Converged For Iregular Positive Energy Solution '//  &
                     'L = '//itoc(angular_momentum)//' in ', count,' Terms'
    END IF
  END DO
  pre = r**(-angular_momentum - 1)
  dl = dl_fun()
  gl = dl*gl*pre*r
  dgl = dl*dgl*pre
!**********************************************************************c
!               continue with rest of wavefunction                     c
!**********************************************************************c
  pre_fac = two * eta_in / c0_sq_fun()
  argm = LOG (two*r) + ql_pl_fun()
  gl = gl + pre_fac * fl * argm
  dgl = dgl + pre_fac * (dfl*argm + fl * r_inv)
  END SUBROUTINE series_expansion_irregular_positive_energy_function
!=======================================================================
!=======================================================================
!deck series_expansion_regular_negative_energy_function
!***begin prologue     series_expansion_regular_negative_energy_function
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            series expansion for negative energy
!***                   regular coulomb function at small rho.
!***description
!***references         NBS handbook
!***routines called
!***end prologue       regular_negative_energy_function
  SUBROUTINE series_expansion_regular_negative_energy_function(a_0)
  IMPLICIT NONE
  REAL*8                     :: a_0( 0 : series_size)
  fl = zero
  dfl = zero
  pre = one
  two_rho =two*r
  count = int_zero
  DO i=0,series_size
     count = count + int_one
     add=a_0(i)*pre
!**********************************************************************c
!            this funny business is done because accidental zeros      c
!            can occur in the a coefficients and at least a            c
!            reasonable result is likely to come out with 25 terms     c
!**********************************************************************c
    IF (ABS(add) < tol.AND.i >= 25) EXIT
    fl = fl + add
    dfl = dfl + i*add
    pre = pre*two_rho
    IF (print_convergence) THEN
        WRITE (iout,*) 'Power Series Converged For Regular Negative Energy Solution '//  &
                     'L = '//itoc(angular_momentum)//' in ', count,' Terms'
    END IF
  END DO
  pre = invsqrt2*(two_rho)**(angular_momentum + 1)
  add = EXP(-r)
  fl = fl*pre*add
  dfl = fl*( -one + (angular_momentum + 1 )*r_inv ) + pre*add*dfl*r_inv
  END SUBROUTINE series_expansion_regular_negative_energy_function
!=======================================================================
!=======================================================================
!deck series_expansion_irregular_negative_energy_function
!***begin prologue     series_expansion_irregular_negative_energy_function
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            series expansion for negative energy
!***                   series_expansion_irregular coulomb function at small rho.
!***description
!***references
!***routines called
!***end prologue       series_expansion_irregular_negative_energy_function
  SUBROUTINE series_expansion_irregular_negative_energy_function(a_0,b_0,c_0)
  IMPLICIT NONE
  REAL*8                                 :: a_0( 0 : series_size)
  REAL*8                                 :: b_0( 0 : series_size)
  REAL*8                                 :: c_0( 0 : series_size)
  INTEGER                                :: i
!**********************************************************************c
!      gl = -1./sqrt(2.) * series * exp(-rho) * (2.*rho)**(-l) /       c
!                    (gamma(eta+l+1)*gamma(eta-l))                     c
!                                                                      c
!      series = sum(n=0 to n=2*l) * c0(n) * (-2.*rho)**n               c
!               - (-2.*rho)**(2*l+1) * { ln(2.*rho) *                  c
!                                        sum(n=0 to n=infinity)        c
!                                            a0(n) * (2.*rho)**n -     c
!                                        sum(n=0 to n=infinity)        c
!                                            b0(n) * (2.*rho)**n }     c
!                                                                      c
!      * Henry has two misprints in his paper. the c0-sum has a        c
!        + 2.*rho and the last sum is added instead of subtracted.     c
!        these errors were discovered by comparing with the NBS        c
!        handbook and using the Henry normalization.                   c
!**********************************************************************c
  two_el = angular_momentum + angular_momentum
  two_el_1 = two_el + int_one
  pre_fac = invsqrt2 / (gamma(eta_in+angular_momentum+1)*gamma(eta_in-angular_momentum))
!**********************************************************************c
!        exponential prefactor and its derivative                      c
!**********************************************************************c
  exp_fac = pre_fac*EXP(-r)
  dxp_fac = -exp_fac
!**********************************************************************c
!             (2.*rho)**(-l) and derivative                            c
!             then combine with exponential factor                     c
!**********************************************************************c
  i_fac = one
  IF (angular_momentum /= 0) THEN
    i_fac = (half*r_inv)**angular_momentum
    exp_fac = i_fac*exp_fac
    dxp_fac = i_fac*dxp_fac-angular_momentum*exp_fac*r_inv
  END IF
!**********************************************************************c
!                prefactors in front of sums                           c
!**********************************************************************c
  two_rho = two*r
  fac = (-two_rho)**two_el_1
  lg_fac = LOG(two_rho)
!**********************************************************************c
!                     do the six sums                                  c
!**********************************************************************c
  c_sum = zero
  dc_sum = zero
  pre = one
  DO i=0,two_el
     c_sum = c_sum+c_0(i)*pre
     dc_sum = dc_sum+i*pre
     pre = - pre * two_rho
  END DO
  dc_sum = dc_sum*r_inv
  pre = one
  a_sum = zero
  da_sum = zero
  DO i=0,series_size
     add = a_0(i)*pre
!**********************************************************************c
!            this funny business is done because accidental zeros      c
!            can occur in the a coefficients and at least a            c
!            reasonable result is likely to come out with 25 terms     c
!**********************************************************************c
     IF (ABS(add) < tol.AND.i >= 25) EXIT
     a_sum = a_sum+add
     da_sum = da_sum+i*add
     pre = pre*two_rho
  END DO
  i_max = i
  da_sum = da_sum*r_inv
  pre = one
  b_sum = zero
  db_sum = zero
  DO i=0,series_size
     add = b_0(i)*pre
     IF (ABS(add) < tol.AND.i >= 25) EXIT
     b_sum = b_sum+add
     db_sum = db_sum+i*add
     pre = pre*two_rho
  END DO
  i_max=MAX(i_max,i)
  IF (print_convergence) THEN
        WRITE (iout,*) 'Power Series Converged For Iregular Negative Energy Solution '//  &
                     'L = '//itoc(angular_momentum)//' in ', i_max,' Terms'
  END IF
  db_sum = db_sum*r_inv
!**********************************************************************c



!        combine all the factors to get the irregular function         c
!                                and                                   c
!                             its   derivative                         c
!**********************************************************************c
  fac_1 = a_sum*lg_fac - b_sum
  fac_2 = a_sum*r_inv + lg_fac*da_sum - db_sum
  fac_3 = c_sum-fac*fac_1
  gl = exp_fac*fac_3
  dgl = exp_fac*(dc_sum-fac*fac_2-two_el_1*fac*fac_1*r_inv) + +dxp_fac*fac_3
  dgl = - dgl
  gl = -gl
  END SUBROUTINE series_expansion_irregular_negative_energy_function
!=======================================================================
!=======================================================================
!deck asymptotic_expansion_positive_energy_function
!***begin prologue     asymptotic_expansion_positive_energy_function
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            asymptotic expansion for regular or irregular
!***                   positive coulomb function at large rho.
!***
!***description        f  = g * cos ( theta ) + f * sin ( theta )
!                       l                  l                   l
!***                   the leading term of f is 1.0
!***                   for the irregular function switch f and g and replace
!***                   the plus by a minus sign.
!***                   g, f and their derivatives are defined in the NBS
!***                   mathematical handbook.
!***references

!***routines called

!***end prologue       asymptotic_expansion_positive_energy_function
  SUBROUTINE asymptotic_expansion_positive_energy_function(a_i,b_i)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                       :: a_i(0:asymptotic_size)
  REAL*8, DIMENSION(:)                       :: b_i(0:asymptotic_size)
!**********************************************************************c
!        calculate the coulomb phase shift and the sin and cos         c
!        of the argument of the asymptotic theta_l                     c
!**********************************************************************c
  sigma_l = c_arg(cgamma(angular_momentum+one+eye*eta_in))
  IF( print_sigma_l ) THEN
     WRITE(iout,1) sigma_l
  END IF
  theta = r - eta_in*LOG(two*r)  - angular_momentum*pi*half 
  theta_l = theta + sigma_l
  s_n(1) = SIN(theta_l)
  s_n(2) = SIN(theta)
  c_n(1) = COS(theta_l)
  c_n(2) = COS(theta)
  s_sigma = sin(sigma_l)
  c_sigma = cos(sigma_l)
!**********************************************************************c
!                do the series expansion to the level tol              c
!**********************************************************************c
  f_old = one
  g_old = zero
  f_old_d = zero
  g_old_d = one - eta_in*r_inv
  f = f_old
  g = g_old
  f_d = f_old_d
  g_d = g_old_d
  DO i=0, asymptotic_size
     f_new = (a_i(i)*f_old-b_i(i)*g_old)*r_inv
     g_new = (a_i(i)*g_old+b_i(i)*f_old)*r_inv
     f_d_new = (a_i(i)*f_old_d-b_i(i)*g_old_d-f_new)*r_inv
     g_d_new = (a_i(i)*g_old_d+b_i(i)*f_old_d-g_new)*r_inv
     IF (ABS(f_new) < tol.AND.ABS(g_new) < tol) EXIT
     f = f+f_new
     g = g+g_new
     f_d = f_d+f_d_new
     g_d = g_d+g_d_new
     f_old = f_new
     g_old = g_new
     f_old_d = f_d_new
     g_old_d = g_d_new
  END DO
  IF ( print_convergence ) THEN
      WRITE(iout,*) 'asymptotic expansion for L = '//itoc(angular_momentum)//' converged in ',i, ' terms'
  END IF
  fl = f*s_n(1)    + g*c_n(1)
  gl = f*c_n(1)    - g*s_n(1)
  dfl = f_d*s_n(1) + g_d*c_n(1)
  dgl = f_d*c_n(1) - g_d*s_n(1)
  phase_factor(1) =  f*s_n(2) + g*c_n(2)
  phase_factor(2) =  f*c_n(2) - g*s_n(2) 
  phase_factor(3) =  f_d*s_n(2) + g_d*c_n(2)
  phase_factor(4) =  f_d*c_n(2) - g_d*s_n(2)
1 FORMAT(/,5X,'coulomb phase shift',1X,e15.8,/)
END SUBROUTINE asymptotic_expansion_positive_energy_function
!=======================================================================
!=======================================================================
!***********************************************************************
!*deck asymptotic_expansion_irregular_negative_energy_function
!***begin prologue     asymptotic_expansion_irregular_negative_energy_function
!***date written       920324   (yymmdd)
!***revision date               (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             %W% %G%
!***purpose            asymntotic expansion for irregular negative
!***                   energy coulomb function at large rho.
!***
!***description        exponentially decaying coulomb function is
!***                   computed using an asymptotic expansion at
!***                   large rho.
!***references

!***routines called

!***end prologue       asymn
  SUBROUTINE asymptotic_expansion_irregular_negative_energy_function(e_0)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                  :: e_0
  REAL*8   pre, add, fac
!**********************************************************************c
!            gl = -1./sqrt(2.) * (series) * exp(-rho)                  c
!                              * (2.rho)**(-eta)                       c
!                                                                      c
!            series = sum ( n =0 to n= infinity ) e0(n)                c
!                             * (-2.*rho)**(-eta)                      c
!                                                                      c
!                * Henry incorrectly uses 2.*rho                       c
!**********************************************************************c
!**********************************************************************c
!                do the series expansion to the level tol              c
!**********************************************************************c
  pre=one
  fac=-half*r_inv
  DO i = 0, asymptotic_size
     add=e_0(i)*pre
     IF (ABS(add) <= tol) EXIT
     gl = gl + add
     dgl = dgl + i*add
     pre=pre*fac
  END DO
  IF ( print_convergence ) THEN
      WRITE(iout,*) 'asymptotic expansion for L = '//itoc(angular_momentum)//' converged in ',i, ' terms'
  END IF
  fac = half*r_inv
  pre = invsqrt2*EXP(-r)
  add = fac**eta_in
  gl = -gl*pre*add
  dgl = -gl*(one+eta_in*r_inv)-pre*add*dgl*r_inv
END SUBROUTINE asymptotic_expansion_irregular_negative_energy_function
!=======================================================================
!=======================================================================
!                                                                      *
!  C O U L F G  -  P A C K A G E   (FROM  A.R. B A R N E T T)          *
!                                                                      *
!  THIS PACKAGE IS USED TO CALCULATE REGULAR AND IRREGULAR             *
!  COULOMB - AND BESSEL - FUNCTIONS                                    *
!                                                                      *
SUBROUTINE coulfg(xx,eta_in,fc,gc,dfc,dgc)
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
!  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA_in (INCLUDING 0), 
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
REAL*8                                   :: eta_in
INTEGER                                  :: mode1
INTEGER                                  :: kfn
INTEGER                                  :: iexp
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
  eta = eta_in
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
END MODULE Coulomb_Functions_Module
