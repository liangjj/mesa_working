MODULE GLpotential_module 

IMPLICIT NONE

! *** CONSTANT PARAMETERS ***
  
	REAL(KIND=8), PARAMETER :: PI = 3.14_8
	REAL(KIND=8), PARAMETER :: HBAR = 1.05459E-34_8
	REAL(KIND=8), PARAMETER :: LINEWIDTH = 2.0_8*PI*9.8E6_8 ! linewidth of excited state
	REAL(KIND=8), PARAMETER :: CLIGHT = 2.998E8_8
	REAL(KIND=8), PARAMETER :: D2DETUNING = -1.5E9_8 ! detuning from D_2 line
	REAL(KIND=8), PARAMETER :: D1DETUNING = 516E9_8 + D2DETUNING ! detuning from D_1 line
	REAL(KIND=8), PARAMETER :: GBEAMWIDTH = 175E-6_8 ! w0 of gaussian beam equal to 1/2 beam diameter 
	REAL(KIND=8), PARAMETER :: GBEAMPOWER = 10E-6_8 ! power in the Gaussian beam
	REAL(KIND=8), PARAMETER :: GLBEAMWIDTH = (120E-6_8)/1.4142_8 ! w0 of LG-beam equal to 1/sqrt(2) times beam diameter
	REAL(KIND=8), PARAMETER :: GLBEAMPOWER = 1.5E-6_8 ! power in LG beam
	REAL(KIND=8), PARAMETER :: GPEAKINTENS = 2.0_8*GBEAMPOWER/(PI*GBEAMWIDTH*GBEAMWIDTH) ! peak intensity in gaussian beam
	REAL(KIND=8), PARAMETER :: OMEGADIFF = 2.0_8*PI*100E3_8 ! frequency difference between the two counter-propagating beams
	REAL(KIND=8), PARAMETER :: WAVEVECTOR = 2.0_8*PI/589.16E-9 ! the wave vector of the light

CONTAINS

	FUNCTION GLpotential(r,phi,z,t)
		! *** INPUT ***
		! r is the radial coordinate, in units of meters
		! phi is the azimuthal coordinate
		! z is the axial coordinate, in units of meters
		! t is the unit of time in units of seconds
		! *** OUTPUT ***
		! GLpotential is the value of the optical potential in units of Hz

		REAL(KIND=8), INTENT(IN) :: r, phi, z, t
		REAL(KIND=8) :: GLpotential, A1, A2, intensity
		A1=sqrt(GPEAKINTENS)*EXP(-r*r/GBEAMWIDTH/GBEAMWIDTH); ! gaussian profile
		A2=sqrt(GLBEAMPOWER)*2.0/(SQRT(PI)*GLBEAMWIDTH*GLBEAMWIDTH)*r*EXP(-r*r/GLBEAMWIDTH/GLBEAMWIDTH); ! LG profile
		intensity=(A1*COS(-WAVEVECTOR*z) + A2*COS(WAVEVECTOR*z + OMEGADIFF*t + phi))**2 &
				 +(A1*SIN(-WAVEVECTOR*z) + A2*SIN(WAVEVECTOR*z + OMEGADIFF*t + phi))**2; 
				 ! the light intensity distribution as a function of space and time 
		GLpotential=intensity*PI*CLIGHT*CLIGHT*LINEWIDTH/2.0_8/(2.0_8*PI*CLIGHT/589E-9)**3 &
				   *(2.0_8/(D2DETUNING*2.0_8*PI) + 1.0_8/(D1DETUNING*2.0_8*PI))/(2.0_8*PI*HBAR); 
				 ! the dipole potential in Hz
	END FUNCTION GLpotential


END MODULE GLpotential_module


