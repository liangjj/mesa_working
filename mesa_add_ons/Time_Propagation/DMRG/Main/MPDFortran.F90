PROGRAM MPDFortran
	USE MPDtools_module
	USE io_module
	USE local_operations_module
	USE observables_module
	USE propagation_module
	IMPLICIT NONE
#ifdef MPI
	INCLUDE 'mpif.h'
	INTEGER :: stat(MPI_STATUS_SIZE)
#endif MPI
	INTEGER :: numNodes, ierr, nodeID, i, j, k, numGridPoints, assignment, master, sender, dummy, imu, iJU, &
	nodeDataIndex, chi, w0, w1, w2
	REAL(KIND=8) :: dJU, dmu, muInput, tunnelInput, U0Input, number, energy, depletion, gap, ncenter
	TYPE(vector) :: tunnel, mu, U0
	INTEGER, ALLOCATABLE :: JUIndex(:), muIndex(:)
        REAL(KIND=8), ALLOCATABLE :: nodeData(:,:)
#ifdef MPI
	CALL mpi_init(ierr)
	CALL mpi_comm_size(MPI_COMM_WORLD, numNodes, ierr)
	CALL mpi_comm_rank(MPI_COMM_WORLD, nodeID, ierr)
#endif MPI
	! *** INITIALIZE PARAMETERS ***
        OPEN(inp,file='DMRG_input',status='old')
        OPEN(iout,file='DMRG_output',status='unknown')
        WRITE(iout,1)
        WRITE(iout,2)
        WRITE(iout,1)
	READ(inp,itpIntParams)
        WRITE(iout,3)  systemSize, maxBoseFilling, ITPswitch, messageSwitchITP, maxITPsteps, chiMin, &
                       chiMax, chiInc, chiConvSwitch, muPoints, JUPoints, JUSwitch
	READ(inp,itpRealParams)
        WRITE(iout,4)  dtITP, numberTolInner, numberTolOuter, muMin, muMax, JUMin, JUMax, trap, lambda
	READ(inp,itpFileParams)
        WRITE(iout,5)  itpDir, itpExt
	READ(inp,evoIntParams)
        WRITE(iout,6) chiEnlargeSwitch, chiEvo, numDataStores, evoSwitch, messageSwitchEvo, xShiftPoints
	READ(inp,evoRealParams)
        WRITE(iout,7)  dtEvo, evolveTime, x0Min, x0Max
	READ(inp,evoFileParams)		
        WRITE(iout,8)  evoInDir, evoInExt, evoOutDir, evoOutExt, iIn, jIn
	localSize=maxBoseFilling+1
	numGridPoints=muPoints*JUPoints
	master=0
	dummy=0
	ALLOCATE(tunnel%v(JUPoints))
	ALLOCATE(mu%v(muPoints))
	ALLOCATE(U0%v(JUPoints))
	ALLOCATE(JUIndex(numGridPoints))
	ALLOCATE(muIndex(numGridPoints))
#ifdef MPI
        ALLOCATE(nodeData(numGridPoints,14))
#else
        ALLOCATE(nodeData(1,14))
#endif MPI
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
           tunnel%v(2:JUPoints) = tunnel%v(1:JUPoints-1) + dJU
           U0%v(2:JUPoints) = 1.0_8
	ELSE
           tunnel%v(1)=1.0_8
           U0%v(1)=JUMin
           tunnel%v(2:JUPoints)=1.0_8
           U0%v(2:JUPoints)=U0%v(1:JUPoints-1) + dJU
	END IF		
	mu%v(1)=muMin
        mu%v(2:muPoints)=mu%v(1:muPoints-1) + dmu
#ifdef MPI
	k=0
	DO i=1,muPoints
           DO j=1,JUPoints
              k=k+1
              JUIndex(k)=j
              muIndex(k)=i
          END DO
       END DO
	PRINT *, nodeID
	IF(nodeID==master) THEN
            ITPParamsName=makeFileName('ITPparams',itpDir,itpExt,nodeID,9)
            progressName=makeFileName('PROGRESS',itpDir,itpExt,nodeID,8)
            CALL OpenFile(3,ITPParamsName)
            WRITE(UNIT=3, FMT=*) systemSize, localSize, muPoints, JUPoints, JUSwitch, &
                                 trap, lambda, maxBoseFilling, maxITPsteps,           &
                                 chiConvSwitch, dtITP, numberTolInner, numberTolOuter
            CLOSE(3)
            DO i=1,(numNodes-1)
               CALL MPI_SEND(i,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,ierr)
            END DO
            OPEN(UNIT=7, FILE=progressName, POSITION="REWIND")
            WRITE(UNIT=7, FMT=*) 'Begin'
            CLOSE(7)
            DO i=numNodes,numGridPoints
               CALL MPI_RECV(dummy,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,          &
                             MPI_COMM_WORLD,stat,ierr)
               sender=stat(MPI_SOURCE)
               CALL MPI_SEND(i,1,MPI_INTEGER,sender,i,MPI_COMM_WORLD,ierr)
               IF(MOD(i,1)==0) THEN
                  OPEN(UNIT=7, FILE=progressName, POSITION="APPEND")
                  WRITE(UNIT=7, FMT=*) i, muIndex(i), JUIndex(i), sender
                  CLOSE(7)
               END IF
            END DO
		DO i=1,(numNodes-1)
                   CALL MPI_SEND(0,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,ierr)
                END DO		
        ELSE
            nodeDataName=makeFileName('nodeData',itpDir,itpExt,nodeID,8)
            assignment=1
            CALL MPI_RECV(assignment,1,MPI_INTEGER,master,MPI_ANY_TAG,                &
                          MPI_COMM_WORLD,stat,ierr)
            nodeDataIndex=0
            DO WHILE (assignment/=0)
               nodeDataIndex=nodeDataIndex+1
               imu=muIndex(assignment)
               iJU=JUIndex(assignment)
               muInput=mu%v(imu)
               tunnelInput=tunnel%v(iJU)
               U0Input=U0%v(iJU)
               CALL PhaseDiagram(muInput, tunnelInput, U0Input, chi, number, energy,  &
                                 depletion, gap, ncenter, w0, w1, w2)
               nodeData(nodeDataIndex,1)=imu
               nodeData(nodeDataIndex,2)=iJU
               nodeData(nodeDataIndex,3)=chi
               nodeData(nodeDataIndex,4)=w0
               nodeData(nodeDataIndex,5)=w1
               nodeData(nodeDataIndex,6)=w2
               nodeData(nodeDataIndex,7)=muInput
               nodeData(nodeDataIndex,8)=tunnelInput
               nodeData(nodeDataIndex,9)=U0Input
               nodeData(nodeDataIndex,10)=number
               nodeData(nodeDataIndex,11)=energy
               nodeData(nodeDataIndex,12)=depletion
               nodeData(nodeDataIndex,13)=gap
               nodeData(nodeDataIndex,14)=ncenter
               CALL MPI_SEND(dummy,1,MPI_INTEGER,master,assignment,MPI_COMM_WORLD,ierr)
               CALL MPI_RECV(assignment,1,MPI_INTEGER,master,MPI_ANY_TAG,             &
                             MPI_COMM_WORLD,stat,ierr)
            END DO
            WRITE(UNIT=7, FMT=*) INT(nodeData(1:nodeDataIndex,1:6)),                  &
                                    (nodeData(1:nodeDataIndex,7:14))            
            CLOSE(7)
        END IF
	DEALLOCATE(JUIndex)
	DEALLOCATE(muIndex)
	CALL mpi_finalize(ierr) 
#else
        DMRG_DataName = makeFileNameSingle('DMRG_Data',itpDir,itpExt,10)
        CALL OpenFile(7,DMRG_DataName)
        WRITE(iout,*) '                             OUTPUT DATA'
        DO i=1, muPoints
           nodeData(1,1)=i
           muInput=mu%v(i)
           DO j=1,JUpoints
              nodeData(1,2)=j
              tunnelInput=tunnel%v(j)
              U0Input=U0%v(j)
              CALL PhaseDiagram(muInput, tunnelInput, U0Input, chi, number, energy,  &
                                depletion, gap, ncenter, w0, w1, w2)
              nodeData(1,3)=chi
              nodeData(1,4)=w0
              nodeData(1,5)=w1
              nodeData(1,6)=w2
              nodeData(1,7)=muInput
              nodeData(1,8)=tunnelInput
              nodeData(1,9)=U0Input
              nodeData(1,10)=number
              nodeData(1,11)=energy
              nodeData(1,12)=depletion
              nodeData(1,13)=gap
              nodeData(1,14)=ncenter
              WRITE(UNIT=7, FMT=*) INT(nodeData(1,1:6)), nodeData(1,7:14)                        
              WRITE(iout,9) INT(nodeData(1,3:6))            
              WRITE(iout,10) nodeData(1,7:14)            
           END DO
        END DO
#endif MPI
        CLOSE(7)
	! *** CLOSE AND DEALLOCATE
	CLOSE(1)
	CALL DestroyFieldOps()
	DEALLOCATE(tunnel%v)
	DEALLOCATE(mu%v)
	DEALLOCATE(U0%v)
	DEALLOCATE(nodeData)
1  FORMAT(/,1x,'************************************************************************')
2  FORMAT(/,20x,'Density Matrix Renormalization Program')
3  FORMAT(/,5x,'system_Size        = ',i5,1x,'max_Bose_Filling = ',i5,1x,'ITP_switch = ',i5,/,5x,           &
              'message_Switch_ITP = ',i5,1x,'max_ITP_steps    = ',i10,                     /,5x,            &
              'chi_Min            = ',i5,1x,'chi_Max          = ',i5,1x,'chi_Inc    = ',i5,/,5x,            &
              'chi_Conv_Switch    = ',i5,1x,'mu_Points        = ',i5,1x,'JU_Points  = ',i5,/,5x,            &
              'JU_Switch          = ',i5)
4  FORMAT(/,5x,'dt_ITP = ',e15.8,1x,'number_Tol_Inner = ',e15.8,1x,'number_Tol_Outer = ',e15.8,/,5x,        &
              'mu_Min = ',e15.8,1x,'mu_Max           = ',e15.8,1x,'JU_Min           = 'e15.8,/,5x,          &
              'JU_Max = ',e15.8,1x,'trap             = ',e15.8,1x,'lambda           = ',e15.8)
5  FORMAT(/,5x,'itp_Dir  = ',a10,1x,'itr_Ext = ',a10)
6  FORMAT(/,5x,'chi_Enlarge_Switch = ',i5,1x,'chi_Evo    = ',i5,                              /,5x,         &
              'num_Data_Stores    = ',i5,1x,'evo_Switch = ',i5,1x,'message_Switch_Evo = ',i5,/,5x,          &
              'x_Shift_Points     = ',i5)
7  FORMAT(/,5x,'dt_Evo = ',e15.8,1x,'evolve_Time = ',e15.8,1x,'x0_Min = ',e15.8,1x,'x0_Max = ',e15.8)
8  FORMAT(/,5x,'evo_In_Dir   = ',a10,1x,'evo_In_Ext  = ',a10,/,5x,                                          &
               'evo_Out_Dir  = ',a10,1x,'evo_Out_Ext = ',a10,/,5x,                                          &
              'i_In = ',i5,1x,'j_In = ',i5)
9  FORMAT(/,5x,'chi = ',i5,1x,'w0 = ',i5,1x,'w1 = ',i5,1x,'w2 = ',i5)
10 FORMAT(/,5x,'mu_Input = ',e15.8,1x,'tunnel_Input = ',e15.8,1x,'U0_Input  = ',e15.8,/,5x,                 &
               'number   = ',e15.8,1x,'energy       = ',e15.8,1x,'depletion = 'e15.8,/,5x,                      &
               'gap      = ',e15.8,1x,'n_center     = ',e15.8)
END PROGRAM MPDFortran









