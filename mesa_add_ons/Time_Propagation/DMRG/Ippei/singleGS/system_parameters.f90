MODULE system_parameters 

IMPLICIT NONE

! *** GLOBAL INPUT PARAMETERS

	INTEGER :: systemSize, maxBoseFilling, ITPswitch, messageSwitchITP, maxITPsteps, chiMin, chiMax, chiInc, chiConvSwitch, &
			   muPoints, JUPoints, JUSwitch
	REAL(KIND=8) :: dtITP, numberTolInner, numberTolOuter, muMin, muMax, JUMin, JUMax, trap, lambda
	CHARACTER(32) :: itpDir, itpExt
	INTEGER :: chiEnlargeSwitch, chiEvo, numDataStores, evoSwitch, messageSwitchEvo, xShiftPoints, iIn, jIn
	REAL(KIND=8) :: dtEvo, evolveTime, x0Min, x0Max
	CHARACTER(32) :: evoInDir, evoInExt, evoOutDir, evoOutExt

! *** GLOBAL DERIVED PARAMETERS
	INTEGER :: localSize
	REAL(KIND=8) :: dtBig


END MODULE system_parameters