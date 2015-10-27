!deck Legendre_Function_Driver
!**begin prologue     DRIVER
!**date written       060711   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source
!**purpose            test code to compute the regular and/or irregular associated Legendre
!***                  functions for integer values of the degree and order and any real argument.
!**description        see CPC paper or type 'yes' at prompt.
!**                   
!**references
!**routines called        
!                       ----                    --------
!
!**end prologue       DRIVER
  PROGRAM Legendre_Function_Driver
  USE accuracy
  USE input_output
  USE Associated_Legendre_Functions
  USE Lentz_Thompson
  IMPLICIT NONE
!                                                                                                  
  TYPE (CF_Legendre)                       :: CFL
  CHARACTER(LEN=8)                         :: itoc
  CHARACTER(LEN=16)                        :: ans
  INTEGER                                  :: i
  INTEGER                                  :: len
  namelist / input_data / title, l_max, m_max, n_points, upper, lower, &
                          Directive, input_values, normalize, Control, &
                          Print_Functions, Print_Wronskian, eps,       &
                          R, recur, test_wron, smallest, biggest
!
!  Get the input and output file numbers which appear in the 
!  Module input_output.f90
!
!  Open the input and output files
!
  write(6,*) '                  Lets begin the calculation'
  write(6,*) '  If you wish (do not wish) some information, type yes(no)'
  write(6,*) '  First time users should type yes'
  read(5,*) ans
  IF ( ans == 'yes' ) THEN
       Call Info
  ELSE IF (ans == 'no' ) THEN
       OPEN(inp,file='Input_Leg',status='old')
       OPEN(iout,file='Output_Leg',status='unknown')
!
!      If you do wish to take the default values, read the data namelist.
!
       READ(inp,*) ans      
       IF (ans == 'input data' ) THEN
           READ(inp,nml=input_data)
       END IF
!
       write(iout,1)
       len=len_trim(title)
       write(iout,2) title(1:len)
       write(iout,1)
       write(iout,3) l_max, m_max, n_points, Directive, Control, &
                     normalize, eps, recur, test_wron
!
!         Compute factorials
!
       ALLOCATE( Factor(0:l_max+m_max) )
       Call Factorials 
!
!
!       If you simply want to read in arbitrary values of the argument set input_values
!       to .true.  Otherwise read in an upper and lower value, a step and a number of
!       points and the code will generate the arguments.  Note that some compilers do not
!       like allocated variables to appear in namelist statements even though the 
!       allocation is done before the variable is read in.  You might have to fix that.
!       The intel compiler is fine with it.
!
        IF (input_values) THEN
            ALLOCATE(x(1:n_points)) 
            READ(inp,*) x(1:n_points)
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
!
!       Here we have the option of using either the Miller algorithm (A) or the continued
!       fraction approach (B).
!
!
!       Print the arguments.
!
        title='Grid'
        Call Print_Matrix(x,iout,title=title)
        IF ( recur == 'Miller' ) THEN
             Leg%D%A%Dir=recur
        ELSE IF ( recur == 'Wronskian' ) THEN   
             Leg%D%B%Dir=recur
        END IF
!
!       Calculate either regular Legendre only, irregular Legendre only or both.
!
        IF ( Directive == 'regular') THEN
             ALLOCATE(Leg%R_LM%F(0:l_max,0:m_max))
        ELSE IF ( Directive == 'irregular') THEN
             ALLOCATE(Leg%I_LM%F(0:l_max,0:m_max))
        ELSE IF( Directive == 'both') THEN
             ALLOCATE(Leg%R_LM%F(0:l_max,0:m_max), Leg%I_LM%F(0:l_max,0:m_max))
        END IF
        IF (Leg%D%B%Dir == 'Wronskian' ) THEN
            ALLOCATE(Leg%R_L%F(0:l_max))
        END IF
!
!       Do the calculation
!
        IF ( Directive == 'regular') THEN
             Call Legendre( R_LM=Leg%R_LM )
             DEALLOCATE(Leg%R_LM%F)
        ELSE IF ( Directive == 'irregular') THEN
             Call Legendre( I_LM=Leg%I_LM )
             DEALLOCATE(Leg%I_LM%F)
        ELSE IF( Directive == 'both') THEN
             Call Legendre( R_LM=Leg%R_LM, I_LM=Leg%I_LM  )
             DEALLOCATE(Leg%R_LM%F, Leg%I_LM%F)
        END IF
        IF (Leg%D%B%Dir == 'Wronskian' ) THEN
            DEALLOCATE(Leg%R_L%F)
        END IF
        DEALLOCATE( Factor )
        CLOSE(inp)
        CLOSE(iout)
  ELSE
        stop
  END IF
1 Format('           **************************************************************' &
         '****************')
2 Format(15x,'Begin Test Calculation = ',a80)
3 Format(/,10x,'Maximum L                             = ',i6,1x,                     &
               'Maximum M                      = ',i6,                               &
         /,10x 'Number of Points                      = ',i6,1x,                     &
               'Type Functions Computed        = ',a10,1x,                           &
         /,10x,'Type Calculation               = ',a24,1x,                           &
               'Normalize Regular Functions on (-1,1) = ',l1,1x,                     &
         /,10x,'Continued Fraction Convergence = ',1pe15.8,1x,                       &
               'Backward Recurrence Method            = ',a24,1x,                    &
        /,10x, 'Test Wronskian                 = ',l1) 
  stop
END PROGRAM  Legendre_Function_Driver
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
  USE Associated_Legendre_Functions
  write(6,*)  'Compute P_lm(x) and Q_lm(x) for all x'                                 
  write(6,*)  'For P_lm and Q_lm with x=(-1,1) upward recursion is stable.'
  write(6,*)  'For P_lm for abs(x) > 1 upward recursion is also stable'
  write(6,*)  'For Q_lm for abs(x) > 1 upward recursion is unstable and it is necessary to use'
  write(6,*)  'downward recursion.  It is straightforward to obtain the starting values for'
  write(6,*)  'forward recursion.  For downward recursion, we have developed two approaches;'
  write(6,*)  'one is a modification of Millers method where the starting value is obtained'
  write(6,*)  'from a continued fraction representation of the ratio of the last two Q_lm.'
  write(6,*)  'The other approach uses the same continued fraction but combines it with the'
  write(6,*)  'known wronskian of P_lm and Q_lm to obtain exact values of the last two Q_lm.'
  write(6,*)  'The first method requires a renormalization of the Q_lm after the downward'
  write(6,*)  'recursion while the second method, first advanced by Segura and Gil'
  write(6,*)  'is self contained.  The Segura and Gil approach has the disadvantage that it'
  write(6,*)  'is necessary to compute P_lm to obtain the wronskian.  '
  write(6,*)  'If P_lm is also needed by the user, this is not a disadvantage.'
  write(6,*)  'Both approaches produce accurate results.'
  write(6,*)  'Note that the continued fraction does converge slowly when the argument'
  write(6,*)  'is near 1.'
  write(6,*)
  write(6,*)   '              Define'
  write(6,*)   '         s_fac = 1 if abs(z) <= 1 s_fac = -1 if abs(z) > 1'
  write(6,*)   '         y = sqrt( s_fac * z )'
  write(6,*)   '                 Initial Values'
  write(6,*)   '         G_LM is either P_LM or Q_LM'    
  write(6,*)   '         P_MM = - s_fac * ( 2*M - 1) *  sqrt ( s_fac * ( 1 - x* x ) )'
  write(6,*)   '                                     * P_(M-1)(M-1)'
  write(6,*)   '         Q_00 = .5 * ln ( abs( ( z + 1) /( z - 1))'
  write(6,*)   '         Q_10 = z * Q_00 - 1'
  write(6,*)   '         Q_01 = - 1 / sqrt ( s_fac * ( 1 - z * z ) )'
  write(6,*)   '         Q_11 = - s_fac * sqrt ( s_fac * ( 1 - z * z ) ) * ( Q_00 + z / ( 1 - z * z ) )'
  write(6,*)   '                 Recurances'
  write(6,*)   
  write(6,*)  '( L + 1 - M) G_(L+1)M = ( 2*L + 1) z G_LM - ( L + M ) G_(L-1)M'
  write(6,*)  ' is the recursion used to step up or down in l depending on whether'
  write(6,*)  ' you begin with the first or last two members of the sequence'
  write(6,*)  ' To step up in M use'
  write(6, *) 'G_LM = 2*M z G_L(M-1)/ y - s_fac * ( L + M ) * ( L - M + 1 ) * G_L(M-2)'
  write(6,*)  
  write(6,*) '                 Input Variables'
  write(6,*) 'l_max = maximum l  m_max = maximum m  n_points = number of points'
  write(6,*) 'upper = highest value of argument lower = lowest value of argument'
  write(6,*) 'Directive = regular,irregular or both input_values = .true.(Generate from data)'
  write(6,*) 'normalize = .true.(Normalize the functions if on cut)' 
  write(6,*) 'eps = convergence criterion for continued fraction'
  write(6,*) 'recur = Miller or Wronskian x = read in values of argument'
  write(6,*) 'Many variables have default values so that the users will obtain'                                 
  write(6,*) 'output even when they are novices on the use of the code'                                 
  write(6,*) '                  Default Values'
  write(6,*) 'l_max = ', l_max,' m_max = ',m_max,' n_points = ',n_points
  write(6,*) 'lower = ',lower,' upper = ',upper
  write(6,*) 'Directive = ',Directive
  write(6,*) 'Control = ',Control 
  write(6,*) 'normalize = ',normalize,' eps = ',eps
  write(6,*) 'recur = ',recur
  write(6,*) '                     Here is a sample input file'
  write(6,*) '&input_data'
  write(6,*) 'title="test_legendre", l_max=50, m_max=5, Print=.true.'
  write(6,*) 'n_points=1, input_values=.true., upper=1.0, Control="compute_functions",'
  write(6,*) 'lower=-1.0, directive="irregular", normalize=.false.,' 
  write(6,*) 'recur="Wronskian", eps=1.d-15, R = 1.4 /'
  write(6,*) '1.0001 /'
  End Subroutine Info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

