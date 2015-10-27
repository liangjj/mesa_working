!deck DRIVER
!**begin prologue     DRIVER
!**date written       060711   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source
!**purpose            DRIVER test code
!**description        
!**                   
!**references
!**routines called        
!                       ----                    --------
!
!**end prologue       DRIVER
  PROGRAM DRIVER
  USE input_output
  USE Associated_Legendre_Functions
  USE Lentz_Thompson
  IMPLICIT NONE
  TYPE (CF_Legendre)                       :: CFL
  REAL(idp)                                :: upper
  REAL(idp)                                :: lower
  REAL(idp)                                :: step
  LOGICAL                                  :: input_values
  CHARACTER(LEN=8)                         :: itoc
  CHARACTER(LEN=16)                        :: Directive
  CHARACTER(LEN=24)                        :: recur
  CHARACTER(LEN=3)                         :: ans
  INTEGER                                  :: i
  namelist / input_data / title, l_max, m_max, n_points, upper, lower, &
                          directive, input_values, normalize, eps, recur, x
!
!  Get the input and output file numbers which appear in the 
!  Module input_output.f90
!
!  Open the input and output files
!
  write(6,*) '                  Lets begin the calculation'
  write(6,*) '  If you wish (do not wish) some information, type yes(no)'
  read(5,*) ans
  IF ( ans == 'yes' ) THEN
       Call Info
  ELSE IF (ans == 'no' ) THEN
       OPEN(inp,file='Input',status='old')
       OPEN(iout,file='Output',status='unknown')
!
!      Read all data except arguments from namelist
!
       READ(inp,nml=input_data)
!
!      If you simply want to read in arbitrary values of the argument set input_values
!      to .true.  Otherwise read in an upper and lower value, a step and a number of
!      points and the code will generate the arguments.  Note that some compilers do not
!      like allocated variables to appear in namelist statements even though the 
!      allocation is done before the variable is read in.  You might have to fix that.
!      The intel compiler is fine with it.
!
     IF (input_values) THEN
         ALLOCATE(x(1:n_points)) 
         READ(inp,nml=input_data)
         upper = zero
         lower = zero
         DO i = 1, n_points
            upper = max(upper,x(i))
            lower = min(lower,x(i))
         END DO
     ELSE
         step = (upper - lower ) / n_points
         n_points = n_points + int_one
         ALLOCATE(x(1:n_points)) 
         x(1) = lower
         DO i = 2, n_points
            x(i) = x(i-int_one) + step
         END DO
     END IF
     write(iout,1)
     write(iout,2) title
     write(iout,1)
     write(iout,3) l_max, m_max, n_points, lower, upper, directive, &
                   normalize, eps, recur
!
!    Print the arguments.
!
     title='Grid'
     Call Print_Matrix(x,iout)
!
!    Here we have the option of using either the Miller algorithm (A) or the continued
!    fraction approach (B).
!
     IF ( recur == 'Miller' ) THEN
          Leg%D%A%Dir=recur
     ELSE IF ( recur == 'Continued_Fraction' ) THEN   
          Leg%D%B%Dir=recur
     END IF
!
!    Calculate either regular Legendre only, irregular Legendre only or both.
!
     IF ( directive == 'regular') THEN
          Call Legendre( R_LM=Leg%R_LM )
     ELSE IF ( directive == 'irregular') THEN
          Call Legendre( I_LM=Leg%I_LM )
     ELSE IF( directive == 'both') THEN
          Call Legendre( R_LM=Leg%R_LM, I_LM=Leg%I_LM  )
     END IF
     CLOSE(inp)
     CLOSE(iout)
  ELSE
     stop
  END IF
!
1 Format('           **************************************************************' &
         '****************')
2 Format(15x,'Begin Test Calculation = ',a80)
3 Format(/,10x,'Maximum L                             = ',i6,1x,                     &
               'Maximum M                      = ',i6,                               &
         /,10x 'Number of Points                      = ',i6,1x,                     &
               'Lower Limit of Argument        = ',e15.8,                            &
         /,10x,'Upper Limit of Argument               = ',e15.8,1x,                  &
               'Type Functions Computed        = ',a10,                              &
         /,10x,'Normalize Regular Functions on (-1,1) = ',l1,1x,                     &
               'Continued Fraction Convergence = ',e15.8,                            &
        /,10x, 'Backward Recurrance Method            = ',a24) 
  stop
END PROGRAM DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck Info
!**begin prologue     Info
!**date written       060711   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source
!**purpose            Info
!**description        
!**                   
!**references
!**routines called        
!                       ----                    --------
!
!**end prologue       Info
  Subroutine Info
  write(6,*)  'Compute P_lm(x) and Q_lm(x) for all x'                                 
  write(6,*)  'For P_lm and Q_lm with x=(-1,1) upward recursion is stable.'
  write(6,*)  'For P_lm for abs(x) > 1 upward recursion is also stable'
  write(6,*)  'For Q_lm for abs(x) > 1 upward recursion is unstable and we use'
  write(6,*)  'downward recursion starting either at a high value of l and then'
  write(6,*)  'renormalizing using the analytically known first function or using'
  write(6,*)  'the continued fraction approach advocated by Segura and Gil.'
  write(6,*)  'This first method, known as Millers algorithm suffers from the'
  write(6,*)  'defect that one has to have a good estimate of the starting value of'
  write(6,*)  'l which is not easy to estimate.  The second approach uses a continued'
  write(6,*)  'fraction, the known wronskian of P_lm and Q_lm and upward recursion'
  write(6,*)  'for P_lm to get the starting values of Q_lm.  Downward recursion in l'
  write(6,*)  'gets the Q_lm for m = (0,1) and then upward recursion in m gets the'
  write(6,*)  'others.  This is accurate and self contained.  Note that the continued'
  write(6,*)  'fraction does converge slowly when the argument is near 1 but the Miller'
  write(6,*)  'algorithm is even worse since one does not have a clue where to start to'
  write(6,*)  'get accurate values.'
  write(6,*)
  write(6,*)  '( L + 1 - M) P_(L+1)M = ( 2*L + 1) z P_LM - ( L + M ) P_(L-1)M'
  write(6,*)  ' is the recursion used to step up or down in l depending on whether'
  write(6,*)  ' you begin with the first or last two members of the sequence'
  write(6,*)  ' To step up in M use'
  write(6, *) 'P_LM = 2*M z P_L(M-1)/ y - s_fac * ( L + M ) * ( L - M + 1 ) * P_L(M-2)'
  write(6,*)   '       s_fac = 1 if abs(z) <= 1 s_fac = -1 if abs(z) > 1'
  write(6,*)   '       y = sqrt( s_fac * z )'
  write(6,*)  'Input Variables'
  write(6,*) 'l_max = maximum l  m_max = maximum m  n_points = number of points'
  write(6,*) 'upper = highest value of argument lower = lowest value of argument'
  write(6,*) 'directive = regular,irregular or both input_values = .true.(Generate from data)'
  write(6,*) 'normalize = .true.(Normalize the functions if on cut)' 
  write(6,*) 'eps = convergence criterion for continued fraction'
  write(6,*) 'recur = Miller or Continued_Fraction x = read in values of argument'
  End Subroutine Info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

