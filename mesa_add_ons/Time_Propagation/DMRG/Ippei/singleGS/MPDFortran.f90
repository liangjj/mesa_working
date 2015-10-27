PROGRAM MPDFortran

	USE system_parameters
	USE MPDtools_module
	USE io_module
	USE local_operations_module
	USE observables_module
	USE propagation_module

	IMPLICIT NONE

!  	INCLUDE 'mpif.h'

	INTEGER :: numNodes, ierr, nodeID, i, j, k, numGridPoints, assignment, master, sender, dummy, imu, iJU, &
			   nodeDataIndex, chi, w0, w1, w2
!	INTEGER :: stat(MPI_STATUS_SIZE)
	REAL(KIND=8) :: dJU, dmu, muInput, tunnelInput, U0Input, number, energy, depletion, gap, ncenter
	TYPE(vector) :: tunnel, mu, U0
	INTEGER, ALLOCATABLE :: JUIndex(:), muIndex(:)
	REAL(KIND=8), ALLOCATABLE :: nodeData(:,:), singleData(:)
	CHARACTER(32) :: ITPParamsName, progressName, nodeDataName, singleDataName
	
	! *** INPUT PARAMETERS ***
	NAMELIST /itpIntParams/ systemSize, maxBoseFilling, ITPswitch, messageSwitchITP, maxITPsteps, chiMin, chiMax, chiInc, &
							chiConvSwitch, muPoints, JUPoints, JUSwitch
	NAMELIST /itpRealParams/ dtITP, numberTolInner, numberTolOuter, muMin, muMax, JUMin, JUMax, trap, lambda
	NAMELIST /itpFileParams/ itpDir, itpExt
	NAMELIST /evoIntParams/ chiEnlargeSwitch, chiEvo, numDataStores, evoSwitch, messageSwitchEvo, xShiftPoints
	NAMELIST /evoRealParams/ dtEvo, evolveTime, x0Min, x0Max
	NAMELIST /evoFileParams/ evoInDir, evoInExt, iIn, jIn, evoOutDir, evoOutExt

!    CALL mpi_init(ierr)
!	CALL mpi_comm_size(MPI_COMM_WORLD, numNodes, ierr)
!   CALL mpi_comm_rank(MPI_COMM_WORLD, nodeID, ierr)

	! *** INITIALIZE PARAMETERS ***
	OPEN(1,file='inputfile')
	READ(1,itpIntParams)
	READ(1,itpRealParams)
	READ(1,itpFileParams)
	READ(1,evoIntParams)
	READ(1,evoRealParams)
	READ(1,evoFileParams)		
	localSize=maxBoseFilling+1
	numGridPoints=muPoints*JUPoints
	master=0
	dummy=0
	ALLOCATE(tunnel%v(JUPoints))
	ALLOCATE(mu%v(muPoints))
	ALLOCATE(U0%v(JUPoints))
	ALLOCATE(JUIndex(numGridPoints))
	ALLOCATE(muIndex(numGridPoints))
	ALLOCATE(singleData(14))
!	ALLOCATE(nodeData(numGridPoints,14))
	CALL CreateFieldOps()
	
	IF(muPoints/=1) THEN
		dmu=(muMax-muMin)/(muPoints-1)
	ELSE
		dmu=0.0_8
	END IF
	IF(JUPoints/=1) THEN
		dJU=(JUMax-JUMin)/(JUPoints-1)
	ELSE
		dJU=0.0_8
	END IF
	IF(JUSwitch==0) THEN
		tunnel%v(1)=JUMin
		U0%v(1)=1.0_8
		DO i=2,JUPoints
			tunnel%v(i)=tunnel%v(i-1)+dJU
			U0%v(i)=1.0_8
		END DO
	ELSE
		tunnel%v(1)=1.0_8
		U0%v(1)=JUMin
		DO i=2,JUPoints
			tunnel%v(i)=1.0_8
			U0%v(i)=U0%v(i-1)+dJU
		END DO
	END IF		

	mu%v(1)=muMin
	DO i=2,muPoints
		mu%v(i)=mu%v(i-1)+dmu
	END DO

	k=0
	DO i=1,muPoints
		DO j=1,JUPoints
			k=k+1
			JUIndex(k)=j
			muIndex(k)=i
		END DO
	END DO
		
!	PRINT *, nodeID
		
!	IF(nodeID==master) THEN

!		ITPParamsName=makeFileName('ITPparams',itpDir,itpExt,nodeID,9)
!		progressName=makeFileName('PROGRESS',itpDir,itpExt,nodeID,8)
!		CALL OpenFile(3,ITPParamsName)
!		WRITE(UNIT=3, FMT=*) systemSize, localSize, muPoints, JUPoints, JUSwitch, trap, lambda, maxBoseFilling, maxITPsteps, &
!							 chiConvSwitch, dtITP, numberTolInner, numberTolOuter
!		CLOSE(3)
!		DO i=1,(numNodes-1)
!			CALL MPI_SEND(i,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,ierr)
!		END DO

!		OPEN(UNIT=5, FILE=progressName, POSITION="REWIND")
!		WRITE(UNIT=5, FMT=*) 'Begin'
!		CLOSE(5)
		
!		DO i=numNodes,numGridPoints
!			CALL MPI_RECV(dummy,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
!			sender=stat(MPI_SOURCE)
!			CALL MPI_SEND(i,1,MPI_INTEGER,sender,i,MPI_COMM_WORLD,ierr)
!			IF(MOD(i,1)==0) THEN
!				OPEN(UNIT=5, FILE=progressName, POSITION="APPEND")
!				WRITE(UNIT=5, FMT=*) i, muIndex(i), JUIndex(i), sender
!				CLOSE(5)
!			END IF
!		END DO
!		DO i=1,(numNodes-1)
!			CALL MPI_SEND(0,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,ierr)
!		END DO		
!	ELSE
!		nodeDataName=makeFileName('nodeData',itpDir,itpExt,nodeID,8)
		singleDataName = makeFileNameSingle('singleData',itpDir,itpExt,10)
!		assignment=1
!		CALL MPI_RECV(assignment,1,MPI_INTEGER,master,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
!		nodeDataIndex=0
!		DO WHILE (assignment/=0)
!			nodeDataIndex=nodeDataIndex+1
!			imu=muIndex(assignment)
!			iJU=JUIndex(assignment)
!			muInput=mu%v(imu)
!			tunnelInput=tunnel%v(iJU)
!			U0Input=U0%v(iJU)
			muInput = muMin
			tunnelInput = JUMin
			U0Input = 1.0
			CALL PhaseDiagram(muInput, tunnelInput, U0Input, chi, number, energy, depletion, gap, ncenter, w0, w1, w2)
!			nodeData(nodeDataIndex,1)=imu
!			nodeData(nodeDataIndex,2)=iJU
!			nodeData(nodeDataIndex,3)=chi
!			nodeData(nodeDataIndex,4)=w0
!			nodeData(nodeDataIndex,5)=w1
!			nodeData(nodeDataIndex,6)=w2
!			nodeData(nodeDataIndex,7)=muInput
!			nodeData(nodeDataIndex,8)=tunnelInput
!			nodeData(nodeDataIndex,9)=U0Input
!			nodeData(nodeDataIndex,10)=number
!			nodeData(nodeDataIndex,11)=energy
!			nodeData(nodeDataIndex,12)=depletion
!			nodeData(nodeDataIndex,13)=gap
!			nodeData(nodeDataIndex,14)=ncenter
!			CALL MPI_SEND(dummy,1,MPI_INTEGER,master,assignment,MPI_COMM_WORLD,ierr)
!			CALL MPI_RECV(assignment,1,MPI_INTEGER,master,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
!		END DO
			singleData(1)=1
			singleData(2)=1
			singleData(3)=chi
			singleData(4)=w0
			singleData(5)=w1
			singleData(6)=w2
			singleData(7)=muInput
			singleData(8)=tunnelInput
			singleData(9)=U0Input
			singleData(10)=number
			singleData(11)=energy
			singleData(12)=depletion
			singleData(13)=gap
			singleData(14)=ncenter
		CALL OpenFile(7,singleDataName)
!		DO i=1,nodeDataIndex
			WRITE(UNIT=7, FMT=*) (INT(singleData(j)), j=1,6), (singleData(j), j=7,14)
!		END DO
		CLOSE(7)
!	END IF


	! *** CLOSE AND DEALLOCATE
	CLOSE(1)
	CALL DestroyFieldOps()
	DEALLOCATE(tunnel%v)
	DEALLOCATE(mu%v)
	DEALLOCATE(U0%v)
	DEALLOCATE(JUIndex)
	DEALLOCATE(muIndex)
	DEALLOCATE(singleData)
!	DEALLOCATE(nodeData)

!	CALL mpi_finalize(ierr) 


END PROGRAM MPDFortran









