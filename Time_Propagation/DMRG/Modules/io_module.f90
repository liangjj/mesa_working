MODULE io_module
	USE system_parameters_module
	USE MPDtools_module
	USE local_operations_module
	IMPLICIT NONE
        INTEGER                                  :: inp=5       
        INTEGER                                  :: iout=6       
        CHARACTER(32) :: ITPParamsName, progressName, nodeDataName, DMRG_DataName
        NAMELIST /itpIntParams/ systemSize, maxBoseFilling, ITPswitch, messageSwitchITP, maxITPsteps, &
                                chiMin, chiMax, chiInc, chiConvSwitch, muPoints, JUPoints, JUSwitch
        NAMELIST /itpRealParams/ dtITP, numberTolInner, numberTolOuter, muMin, muMax, JUMin, JUMax, trap, lambda
        NAMELIST /itpFileParams/ itpDir, itpExt
        NAMELIST /evoIntParams/ chiEnlargeSwitch, chiEvo, numDataStores, evoSwitch, messageSwitchEvo, xShiftPoints
        NAMELIST /evoRealParams/ dtEvo, evolveTime, x0Min, x0Max
        NAMELIST /evoFileParams/ evoInDir, evoInExt, iIn, jIn, evoOutDir, evoOutExt

                              CONTAINS


	FUNCTION makeFileName(basename,itpDir,itpExt,node,baselen)
		INTEGER, INTENT(IN) :: baselen, node
		CHARACTER(baselen), INTENT(IN) :: basename
		CHARACTER(32), INTENT(IN) :: itpDir, itpExt
		CHARACTER(32) :: makeFileName
		CHARACTER(6) :: filesuffix
		CHARACTER*1 :: x
		filesuffix(1:)='_'
		WRITE(x,'(I1)')FLOOR(node/10.)
		filesuffix(2:)=x
		WRITE(x,'(I1)')MOD(node,10)
		filesuffix(3:)=x
		makeFileName=itpDir
		makeFileName(LEN_TRIM(makeFileName)+1:)=basename
		makeFileName(LEN_TRIM(makeFileName)+1:)=itpExt
		makeFileName(LEN_TRIM(makeFileName)+1:)=filesuffix
	END FUNCTION

        FUNCTION makeFileNameSingle(basename,itpDir,itpExt,baselen)
                INTEGER, INTENT(IN) :: baselen
                CHARACTER(baselen), INTENT(IN) :: basename
                CHARACTER(32), INTENT(IN) :: itpDir, itpExt
                CHARACTER(32) :: makeFileNameSingle
                makeFileNameSingle=itpDir
                makeFileNameSingle(LEN_TRIM(makeFileNameSingle)+1:)=basename
                makeFileNameSingle(LEN_TRIM(makeFileNameSingle)+1:)=itpExt
        END FUNCTION

	SUBROUTINE OpenFile(FileID, filename)    
    	INTEGER, INTENT(IN) :: FileID
		CHARACTER(32), INTENT(IN) :: filename
    	INTEGER :: Status
    	OPEN(UNIT=FileID, FILE=filename, RECL=32384, IOSTAT=Status)
    	IF(Status>0) STOP "*** Cannot open the file ***"
	END SUBROUTINE OpenFile

	SUBROUTINE RecordLambdas(fileid, Lambdas)
		TYPE(vector), POINTER :: Lambdas(:)	
	    INTEGER, INTENT(IN) :: fileid
	    INTEGER :: i, j, chi
		chi = SIZE(Lambdas(2)%v)
	    DO i=1,(systemSize+1)
	    	WRITE(UNIT=fileid, FMT='(100e25.16)') Lambdas(i)%v(1:chi)
	    END DO
	END SUBROUTINE RecordLambdas

	SUBROUTINE RecordGammas(fileid, Gammas)
		TYPE(tensor), POINTER :: Gammas(:)	
		INTEGER, INTENT(IN) :: fileid
		INTEGER :: i,j,k,l,chi
		chi=SIZE(Gammas(1)%t,3)
		DO j=1,chi
			WRITE(UNIT=fileid, FMT='(50e25.16)') (REAL(Gammas(1)%t(1,k,j)), k=1,localSize)
		END DO
		DO l=2,(systemSize-1)
			DO i=1,chi
				DO j=1,chi
					WRITE(UNIT=fileid, FMT='(50e25.16)') (REAL(Gammas(l)%t(i,k,j)), k=1,localSize)
				END DO
			END DO
		END DO
		DO i=1,chi
			WRITE(UNIT=fileid, FMT='(50e25.16)') (REAL(Gammas(systemSize)%t(i,k,1)), k=1,localSize)
		END DO
		DO j=1,chi
			WRITE(UNIT=fileid, FMT='(50e25.16)') (AIMAG(Gammas(1)%t(1,k,j)), k=1,localSize)
		END DO
		DO l=2,(systemSize-1)
			DO i=1,chi
				DO j=1,chi
					WRITE(UNIT=fileid, FMT='(50e25.16)') (AIMAG(Gammas(l)%t(i,k,j)), k=1,localSize)
				END DO
			END DO
		END DO
		DO i=1,chi
			WRITE(UNIT=fileid, FMT='(50e25.16)') (AIMAG(Gammas(systemSize)%t(i,k,1)), k=1,localSize)
		END DO
	END SUBROUTINE RecordGammas
	
	SUBROUTINE ReadGammaLambda(itpParamsFileID, gammaFileID, lambdaFileID, Gammas, Lambdas, chiOut, trapOut, tunnelOut, muOut)
    	INTEGER, INTENT(IN) :: itpParamsFileID, gammaFileID, lambdaFileID
		INTEGER, INTENT(OUT) :: chiOut
		REAL(KIND=8), INTENT(OUT) :: trapOut, tunnelOut, muOut
		TYPE(tensor), POINTER :: Gammas(:)
		TYPE(vector), POINTER :: Lambdas(:)
		TYPE(tensor), POINTER :: GammasIn(:)
		TYPE(vector), POINTER :: LambdasIn(:)
		COMPLEX(KIND=8) :: eye
		REAL(KIND=8), ALLOCATABLE :: gammatemp(:)
	    INTEGER :: i,j,k,l,n,d,chi
		eye = CMPLX(0.0,1.0,KIND=8)
		READ(UNIT=itpParamsFileID, FMT=*) n, d, chi, trapOut, tunnelOut, muOut
		IF(n/=systemSize) THEN
			PRINT *, 'systemSize parameter conflicts with input data'
		END IF
		IF(d/=localSize) THEN
			PRINT *, 'localSize parameter conflicts with input data'
		END IF
		ALLOCATE(gammatemp(localSize))
		CALL AllocateGamLam(GammasIn, LambdasIn, chi)
	    DO i=1,(n+1)
	    	READ(UNIT=lambdaFileID, FMT='(100e25.16)') (LambdasIn(i)%v(j), j=1,chi)
	    END DO
		DO j=1,chi
			READ(UNIT=gammaFileID, FMT='(50e25.16)') (gammatemp(k), k=1,localSize)
			DO k=1,localSize
				GammasIn(1)%t(1,k,j)=gammatemp(k)
			END DO
		END DO
		DO l=2,(systemSize-1)
			DO i=1,chi
				DO j=1,chi
					READ(UNIT=gammaFileID, FMT='(50e25.16)') (gammatemp(k), k=1,localSize)
					DO k=1,localSize
						GammasIn(l)%t(i,k,j)=gammatemp(k)
					END DO
				END DO
			END DO
		END DO
		DO i=1,chi
			READ(UNIT=gammaFileID, FMT='(50e25.16)') (gammatemp(k), k=1,localSize)
			DO k=1,localSize
				GammasIn(systemSize)%t(i,k,1)=gammatemp(k)
			END DO
		END DO
		DO j=1,chi
			READ(UNIT=gammaFileID, FMT='(50e25.16)') (gammatemp(k), k=1,localSize)
			DO k=1,localSize
				GammasIn(1)%t(1,k,j)=GammasIn(1)%t(1,k,j)+eye*gammatemp(k)
			END DO
		END DO
		DO l=2,(systemSize-1)
			DO i=1,chi
				DO j=1,chi
					READ(UNIT=gammaFileID, FMT='(50e25.16)') (gammatemp(k), k=1,localSize)
					DO k=1,localSize
						GammasIn(l)%t(i,k,j)=GammasIn(l)%t(i,k,j)+eye*gammatemp(k)
					END DO
				END DO
			END DO
		END DO
		DO i=1,chi
			READ(UNIT=gammaFileID, FMT='(50e25.16)') (gammatemp(k), k=1,localSize)
			DO k=1,localSize
				GammasIn(systemSize)%t(i,k,1)=GammasIn(systemSize)%t(i,k,1)+eye*gammatemp(k)
			END DO
		END DO
		IF(chiEnlargeSwitch==1) THEN
			CALL AllocateGamLam(Gammas, Lambdas, chiEvo)
			chiOut=chiEvo
		ELSE
			CALL AllocateGamLam(Gammas, Lambdas, chi)
			chiOut=chi
		END IF
		CALL CopyGamLam(Gammas,Lambdas,GammasIn,LambdasIn)
		CALL SVDInit(chiOut)
		DO k=1,(n-1)
			CALL CanonicalForm(Lambdas(k)%v,Gammas(k)%t,Lambdas(k+1)%v,Gammas(k+1)%t,Lambdas(k+2)%v)
		END DO
		DO k=(n-1),1,(-1)
			CALL CanonicalForm(Lambdas(k)%v,Gammas(k)%t,Lambdas(k+1)%v,Gammas(k+1)%t,Lambdas(k+2)%v)
		END DO
		CALL DeallocateGamLam(GammasIn, LambdasIn)
		DEALLOCATE(gammatemp)
		CALL SVDFinish()
	END SUBROUTINE ReadGammaLambda

	SUBROUTINE RecordRho(fileid, rho1)
		INTEGER, INTENT(IN) :: fileid
		TYPE(matrix), POINTER :: rho1(:)
		INTEGER :: i,j,k
	    DO k=1,systemSize
	    	DO i=1,localSize
		    	WRITE(UNIT=fileid, FMT='(5001(e16.8))') (REAL(rho1(k)%m(i,j)), j=1,localSize)
		    END DO
		END DO
		DO k=1,systemSize
	    	DO i=1,localSize
		    	WRITE(UNIT=fileid, FMT='(5001(e16.8))') (AIMAG(rho1(k)%m(i,j)), j=1,localSize)
		    END DO
		END DO		
	END SUBROUTINE RecordRho

	SUBROUTINE RecordTwoSiteObservable(fileid, observable)
		INTEGER, INTENT(IN) :: fileid
		COMPLEX(KIND=8), INTENT(IN) :: observable(:,:)
		INTEGER :: i,j
	    DO i=1,systemSize
		    WRITE(UNIT=fileid, FMT='(5001(e16.8))') (REAL(observable(i,j)), j=1,systemSize)
		END DO
		DO i=1,systemSize
	    	WRITE(UNIT=fileid, FMT='(5001(e16.8))') (AIMAG(observable(i,j)), j=1,systemSize)
		END DO		
	END SUBROUTINE RecordTwoSiteObservable

	
END MODULE io_module
