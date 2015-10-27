MODULE propagation_module

	USE system_parameters_module
	USE MPDtools_module
	USE local_operations_module
	USE observables_module
	USE io_module

	IMPLICIT NONE

CONTAINS

	SUBROUTINE TrotterStep2ndOrder(Udt,Gammas,Lambdas)
		TYPE(matrix), POINTER :: Udt(:)
		TYPE(tensor), POINTER :: Gammas(:)
		TYPE(vector), POINTER :: Lambdas(:)
		INTEGER :: i
		DO i=1,(systemSize-1),2
			CALL TwoSiteOp(i,Udt(i)%m,Gammas,Lambdas)
		END DO
		DO i=2,(systemSize-1),2
			CALL TwoSiteOp(i,Udt(i)%m,Gammas,Lambdas)
		END DO
		DO i=1,(systemSize-1),2
			CALL TwoSiteOp(i,Udt(i)%m,Gammas,Lambdas)
		END DO
	END SUBROUTINE TrotterStep2ndOrder

	SUBROUTINE ImagTimeProp(H, GammasOuter, LambdasOuter, chiIn)
		INTEGER, INTENT(INOUT) :: chiIn
		TYPE(matrix), POINTER :: H(:)
		TYPE(tensor), POINTER :: GammasOuter(:), GammasInner(:)
		TYPE(vector), POINTER :: LambdasOuter(:), LambdasInner(:)
		TYPE(matrix), POINTER :: Uitp(:)
		COMPLEX(KIND=8) :: eye, dtodd, dteven
		INTEGER :: i, imax, j, k, alpha, beta, gamma, chi
		REAL(KIND=8) :: numberInner1, numberInner2, numberOuter1, numberOuter2, testInnerTol, testOuterTol, energy
		CALL AllocateOps(Uitp,systemSize-1,localSize*localSize)
		eye=CMPLX(0.0,1.0,KIND=8)
		dtodd=-eye*dtITP/2.0
		dteven=-eye*dtITP
		CALL ConstructPropagators(H, Uitp, dtodd, dteven)
		IF(chiConvSwitch==1) THEN
			imax=chiMax
		ELSE
			imax=chiIn
		END IF
		numberOuter1=0.0_8
		DO i=chiIn,imax,chiInc
			chi=i
			IF(messageSwitchITP==1) THEN
				PRINT *, 'chi', chi
			END IF
			CALL SVDInit(chi)
			CALL AllocateGamLam(GammasInner, LambdasInner, chi)
			CALL CopyGamLam(GammasInner, LambdasInner, GammasOuter, LambdasOuter)
			DO k=1,(systemSize-1)
				CALL CanonicalForm(LambdasInner(k)%v,GammasInner(k)%t,LambdasInner(k+1)%v,GammasInner(k+1)%t,LambdasInner(k+2)%v)
			END DO
			DO k=(systemSize-1),1,(-1)
				CALL CanonicalForm(LambdasInner(k)%v,GammasInner(k)%t,LambdasInner(k+1)%v,GammasInner(k+1)%t,LambdasInner(k+2)%v)
			END DO
			numberInner1=0.0_8
			DO j=1,maxITPSteps
				IF(MOD(j,50)==0) THEN
					CALL TotalNumber(numberInner2, GammasInner, LambdasInner)
					testInnerTol = DABS((numberInner1-numberInner2)/numberInner2)
					IF(messageSwitchITP==1) THEN
						PRINT *, 'ITP step j', j, 'central lambda', LambdasInner(CEILING(systemSize/2.0))%v
						PRINT *, 'number', numberInner2
						PRINT *, 'numberTolInner', numberTolInner, 'testInnerTol', testInnerTol
					END IF
					IF(testInnerTol<numberTolInner) EXIT
					numberInner1 = numberInner2
				END IF
				CALL TrotterStep2ndOrder(Uitp, GammasInner, LambdasInner)
				DO k=1,(systemSize-1)
					CALL CanonicalForm(LambdasInner(k)%v,GammasInner(k)%t,LambdasInner(k+1)%v,GammasInner(k+1)%t,LambdasInner(k+2)%v)
				END DO
				DO k=(systemSize-1),1,(-1)
					CALL CanonicalForm(LambdasInner(k)%v,GammasInner(k)%t,LambdasInner(k+1)%v,GammasInner(k+1)%t,LambdasInner(k+2)%v)
				END DO
			END DO
			IF((j==maxITPSteps).AND.(messageSwitchITP==1)) THEN
				PRINT *, 'WARNING: j=maxITPSteps in ImagTimeProp'
			END IF
			CALL DeallocateGamLam(GammasOuter, LambdasOuter)
			CALL AllocateGamLam(GammasOuter, LambdasOuter, chi)
			CALL CopyGamLam(GammasOuter, LambdasOuter, GammasInner, LambdasInner)		
			CALL TotalNumber(numberOuter2, GammasOuter, LambdasOuter)
			testOuterTol = DABS((numberOuter1-numberOuter2)/numberOuter2)
			IF(messageSwitchITP==1) THEN
				PRINT *, 'numberTolOuter', numberTolOuter, 'testOuterTol', testOuterTol
			END IF
			CALL SVDFinish()
			CALL DeallocateGamLam(GammasInner, LambdasInner)		
			IF((testOuterTol<numberTolOuter).AND.(chi>4)) EXIT
			numberOuter1 = numberOuter2
		END DO
		chiIn=chi
		IF(messageSwitchITP==1) THEN
			PRINT *, 'norm', NormMPD(GammasOuter,LambdasOuter)
			CALL TotalEnergy(energy, H, GammasOuter, LambdasOuter)
			PRINT *, 'energy', energy
		END IF		
		CALL DeallocateOps(Uitp,systemSize-1)
	END SUBROUTINE ImagTimeProp

	SUBROUTINE DomainFind(rho, w0, w1, w2)
		INTEGER, INTENT(OUT) :: w0, w1, w2
		COMPLEX(KIND=8), INTENT(IN) :: rho(:,:)
		REAL(KIND=8) :: densi
		INTEGER :: i
		w0=0
		w1=0
		w2=0
		DO i=1,systemSize
			densi=REAL(rho(i,i))
			IF(densi>0.005) THEN
				w0=w0+1
			END IF
			IF(densi>1.005) THEN
				w1=w1+1
			END IF
			IF(densi>2.005) THEN
				w2=w2+1
			END IF
		END	DO
	END SUBROUTINE DomainFind

	SUBROUTINE PhaseDiagram(mu, tunnel, U0, chiOut, numberOut, energyOut, depletionOut, gapOut, ncenter, w0, w1, w2)
		REAL(KIND=8), INTENT(IN) :: mu, tunnel, U0
		INTEGER, INTENT(OUT) :: chiOut, w0, w1, w2
		REAL(KIND=8), INTENT(OUT) :: numberOut, energyOut, depletionOut, gapOut, ncenter
		TYPE(tensor), POINTER :: Gammas(:)
		TYPE(vector), POINTER :: Lambdas(:)
		TYPE(matrix), POINTER :: H(:)
		TYPE(matrix) :: oneBodyRho, densDens
		INTEGER :: chi, icenter
		icenter=CEILING(systemSize/2.0_8)
		
		CALL AllocateOps(H,systemSize-1,localSize*localSize)
		ALLOCATE(oneBodyRho%m(systemSize,systemSize))
		ALLOCATE(densDens%m(systemSize,systemSize))

		chiOut=chiMin
		CALL AllocateGamLam(Gammas, Lambdas, chiOut)
		CALL AllStates(Gammas, Lambdas)
		CALL HamiltonianBoseHubbard(H,trap,tunnel,U0,mu,0.0_8,lambda)
		CALL ImagTimeProp(H, Gammas, Lambdas, chiOut)
		CALL TwoSiteExpVal(oneBodyRho%m, 1, n_op%m, adagxa_op%m, Gammas, Lambdas)
		CALL TwoSiteExpVal(densDens%m, 1, nsq_op%m, nxn_op%m, Gammas, Lambdas)
		ncenter=REAL(oneBodyRho%m(icenter,icenter))
		CALL DomainFind(oneBodyRho%m,w0,w1,w2)
		CALL TotalNumber(numberOut, Gammas, Lambdas)
		CALL TotalEnergy(energyOut, H, Gammas, Lambdas)
		CALL Egap(gapOut, oneBodyRho%m, densDens%m, tunnel, numberOut)
		CALL Qdepletion(depletionOut, oneBodyRho%m, numberOut)

		CALL DeallocateGamLam(Gammas,Lambdas)
		CALL DeallocateOps(H,systemSize-1)
		DEALLOCATE(oneBodyRho%m)
		DEALLOCATE(densDens%m)
	END SUBROUTINE PhaseDiagram

END MODULE propagation_module



















