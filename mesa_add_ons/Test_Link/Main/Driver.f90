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
  TYPE (CF_PI)                             :: CFPI
  TYPE (CF_e)                              :: CFe
  REAL(idp)                                :: fpkey
  REAL(idp)                                :: upper
  REAL(idp)                                :: lower
  REAL(idp)                                :: step
  LOGICAL                                  :: dollar
  LOGICAL                                  :: logkey
  INTEGER                                  :: intkey
  LOGICAL                                  :: input_values
  INTEGER                                  :: len
  INTEGER                                  :: lenth
  CHARACTER(LEN=80)                        :: chrkey
  CHARACTER(LEN=80)                        :: cpass
  CHARACTER(LEN=3200)                      :: card
  CHARACTER(LEN=8)                         :: itoc
  CHARACTER(LEN=16)                        :: Directive
  CHARACTER(LEN=1)                         :: recur
  INTEGER                                  :: input
  INTEGER                                  :: output
  INTEGER                                  :: i
  namelist / input_data / title, l_max, m_max, n_points, upper, lower, directive,  &
                          input_values, Method, Derivative, normalize, eps, recur, x
  COMMON /io/ input, output
!
!  Get the input and output file numbers which appear in the 
!  Module input_output.f90 and put them into input and output which appears
!  in the ONLY common block in the code.  This is needed to pass into library
!  routines and for no other purpose.
!
  input = inp
  output = iout
!
!  Open the input and output files
!
  OPEN(input,file='Input',status='old')
  OPEN(output,file='Output',status='unknown')
!
  READ(inp,nml=input_data)
!  IF ( dollar('$data',card,cpass,inp) )THEN
!       title = chrkey(card,'title','No Title',' ')
!       l_max = intkey(card,'l_max',10,' ')
!       m_max = intkey(card,'m_max',l_max,' ')
!       n_points = intkey(card,'number_of_points',10,' ')
!       upper = fpkey(card,'upper_limit',1.d0,' ')
!       lower = fpkey(card,'lower_limit',-1.d0,' ')
!       directive=chrkey(card,'type','regular',' ')
!       normalize=logkey(card,'normalize',.false.,' ')
!       input_values=logkey(card,'input_values',.false.,' ')
!       Method=chrkey(card,'Method','F90',' ')
!       Derivative = logkey(card,'derivative',.false.,' ')
!       eps=fpkey(card,'precision',epsilon(1.d0),' ')
!       recur=chrkey(card,'Downward_Recurance','A',' ')
!  END IF
  IF (input_values) THEN
      ALLOCATE(x(1:n_points)) 
!      Call fparr(card,'points',x,n_points,' ')
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
  write(iout,3) l_max, m_max, n_points, lower, upper, directive, normalize, Method
  title='Grid'
  Call Print_Matrix(x,iout)
  IF ( Method /= 'Continued_Fraction') THEN
       IF ( recur == 'A' ) THEN
           Leg%D%A%Dir=recur
       ELSE IF ( recur == 'B' ) THEN   
           Leg%D%B%Dir=recur
       END IF
       IF ( directive == 'regular') THEN
            Call Legendre( R_LM=Leg%R_LM )
       ELSE IF ( directive == 'irregular') THEN
            Call Legendre( I_LM=Leg%I_LM )
       ELSE IF( directive == 'both') THEN
            Call Legendre( R_LM=Leg%R_LM, I_LM=Leg%I_LM  )
       END IF
  ELSE
       Call Continued_Fractions(CFe)
       Call Continued_Fractions(CFPI)
  END IF
  CLOSE(input)
  CLOSE(output)
!
1 Format('           **************************************************************' &
         '****************')
2 Format(15x,'Begin Test Calculation = ',a80)
3 Format(/,10x,'Maximum L               = ',i3,1x,'Maximum M                  = ',i3,1x,    &
         /,10x 'Number of Points        = ',i6,1x,'Lower Limit of Argument    = ',e15.8,    &
         /,10x,'Upper Limit of Argument = ',e15.8,1x,'Type                    = ',a10,      &
         /,10x,'Normalize               = ',l1,1x,'Method = ',a32 ) 
  stop
END PROGRAM DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
