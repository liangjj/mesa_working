!deck Legendre_Function_Driver
!**begin prologue     DRIVER
!**date written       060711   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source
!**purpose            test code to compute the regular associated Legendre
!***                  functions for integer values of the degree and order and any real argument.
!**description        see CPC paper or type 'yes' at prompt.
!**                   
!**references
!**routines called        
!                       ----                    --------
!
!**end prologue       DRIVER
  PROGRAM Legendre_Function_Driver
  USE Regular_Associated_Legendre_Functions
  IMPLICIT NONE
!                                                                                                  
  TYPE (Reg_LM)                            :: R_LM
  CHARACTER(LEN=8)                         :: itoc
  CHARACTER(LEN=16)                        :: ans
  INTEGER                                  :: i
  INTEGER                                  :: len
  namelist / input_data / title, l_max, m_max, n_points, upper, lower,       &
                          input_values, normalize, Print_Functions
                          
!
!  Get the input and output file numbers which appear in the 
!  Module input_output.f90
!
!  Open the input and output files
!
  OPEN(inp,file='Input_Leg',status='old')
  OPEN(iout,file='Output_Leg',status='unknown')
!
!      If you do wish to take the default values, read the data namelist.
!
  READ(inp,nml=input_data)
!
  write(iout,1)
  len=len_trim(title)
  write(iout,2) title(1:len)
  write(iout,1)
  write(iout,3) l_max, m_max, n_points, normalize
!
!         Compute factorials
!
  ALLOCATE( Factorial(0:l_max+m_max) )
  Call N_Factorial 
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
!
!       Calculate either regular Legendre only, irregular Legendre only or both.
!
  ALLOCATE(Leg%R_LM%F(0:l_max,0:m_max))
!
!       Do the calculation
!
  Call Legendre( R_LM )
  DEALLOCATE(Leg%R_LM%F, Factorial )
  CLOSE(inp)
  CLOSE(iout)
  stop
1 Format('           **************************************************************'        &
         '****************')
2 Format(15x,'Begin Test Calculation = ',a80)
3 Format(/,10x,'Maximum L                             = ',i6,1x,                            &
               'Maximum M                             = ',i6,                               &
         /,10x 'Number of Points                      = ',i6,1x,                            &
               'Normalize Regular Functions on (-1,1) = ',l1 )
  stop
END PROGRAM  Legendre_Function_Driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

