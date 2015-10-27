MODULE observables_module

	USE system_parameters_module
	USE MPDtools_module
	USE local_operations_module, ONLY: OneSiteOp

	IMPLICIT NONE

CONTAINS

	FUNCTION NormMPD(Gammas,Lambdas)
		TYPE(tensor), POINTER :: Gammas(:)
		TYPE(vector), POINTER :: Lambdas(:)
		REAL(KIND=8) :: NormMPD
		COMPLEX(KIND=8) :: temp
		COMPLEX(KIND=8) :: normKernal(SIZE(Gammas(1)%t,3),SIZE(Gammas(1)%t,3))
		COMPLEX(KIND=8) :: GammaTemp(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3))
		INTEGER :: alpha,beta,gamma,eta,i,n,chi
		chi=SIZE(Gammas(1)%t,3)
		DO alpha=1,chi
			DO beta=1,chi
				normKernal(alpha,beta)=CMPLX(0.0,KIND=8)
				DO i=1,localSize
					normKernal(alpha,beta)=normKernal(alpha,beta)+Lambdas(2)%v(alpha)*CONJG(Gammas(1)%t(1,i,alpha)) &
											*Gammas(1)%t(1,i,beta)*Lambdas(2)%v(beta)
				END DO
			END DO
		END DO
		DO n=2,(systemSize-1)
			DO gamma=1,chi
				DO i=1,localSize
					DO beta=1,chi
						GammaTemp(gamma,i,beta)=CMPLX(0.0,KIND=8)
						DO eta=1,chi
							GammaTemp(gamma,i,beta) = GammaTemp(gamma,i,beta) &
												    + normKernal(gamma,eta)*Gammas(n)%t(eta,i,beta)
						END DO
					END DO
				END DO
			END DO
			DO alpha=1,chi
				DO beta=1,chi
					normKernal(alpha,beta)=CMPLX(0.0,KIND=8)
					DO gamma=1,chi
						DO i=1,localSize
							normKernal(alpha,beta) = normKernal(alpha,beta) &
												   + Lambdas(n+1)%v(alpha)*CONJG(Gammas(n)%t(gamma,i,alpha)) &
												   * GammaTemp(gamma,i,beta)*Lambdas(n+1)%v(beta)
						END DO
					END DO
				END DO
			END DO
		END DO
		temp=CMPLX(0.0,KIND=8);
		DO alpha=1,chi
			DO beta=1,chi
				DO i=1,localSize
					temp = temp+CONJG(Gammas(systemSize)%t(alpha,i,1))*normKernal(alpha,beta) &
							*Gammas(systemSize)%t(beta,i,1)
				END DO
			END DO
		END DO
		NormMPD=DSQRT(REAL(temp,KIND=8))
	END FUNCTION NormMPD
	
	SUBROUTINE FormSingleSiteRho(rho1, Lambda0, Gamma1, Lambda1)
		COMPLEX(KIND=8) :: rho1(:,:)
		REAL(KIND=8) :: Lambda0(:), Lambda1(:)
		COMPLEX(KIND=8) :: Gamma1(:,:,:)
		INTEGER :: chi0, chi1, i, j, alpha, beta
		chi0 = SIZE(Gamma1,1)
		chi1 = SIZE(Gamma1,3)
		DO i=1,localSize
			DO j=1,localSize
				rho1(i,j)=CMPLX(0.0,KIND=8)
				DO alpha=1,chi0
					DO beta=1,chi1
						rho1(i,j) = rho1(i,j)+Lambda0(alpha)*Gamma1(alpha,i,beta)*Lambda1(beta) &
										   *Lambda1(beta)*CONJG(Gamma1(alpha,j,beta))*Lambda0(alpha)
					END DO
				END DO
			END DO
		END DO
	END SUBROUTINE FormSingleSiteRho
		
	SUBROUTINE SingleSiteDensityMatrix(rho,Gammas,Lambdas)
		TYPE(matrix), INTENT(OUT) :: rho(:)
		TYPE(tensor), POINTER :: Gammas(:)
		TYPE(vector), POINTER :: Lambdas(:)
		INTEGER :: i
		DO i=1,systemSize
			CALL FormSingleSiteRho(rho(i)%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
		END DO
	END SUBROUTINE SingleSiteDensityMatrix


	SUBROUTINE ThetaKernal(Theta,Lambda1,Gamma,Lambda2)
		COMPLEX(KIND=8), INTENT(OUT) :: Theta(:,:,:,:)
		COMPLEX(KIND=8), INTENT(IN) :: Gamma(:,:,:)
		REAL(KIND=8), INTENT(IN) :: Lambda1(:), Lambda2(:)
		INTEGER :: alpha, beta, eta, i, j
		DO alpha=1,SIZE(Gamma,1)
			DO beta=1,SIZE(Gamma,1)
				DO i=1,localSize
					DO j=1,localSize
						Theta(alpha,i,j,beta)=CMPLX(0.0,KIND=8)
						DO eta=1,SIZE(Gamma,3)
							Theta(alpha,i,j,beta)=Theta(alpha,i,j,beta) &
										 		 +Lambda1(alpha)*Gamma(alpha,i,eta)*Lambda2(eta) &
										 		 *Lambda2(eta)*CONJG(Gamma(beta,j,eta))*Lambda1(beta)
						END DO
					END DO
				END DO
			END DO
		END DO
	END SUBROUTINE ThetaKernal
	
	SUBROUTINE ThetaNext(Theta, Gamma, GammaP, Lambda)
		COMPLEX(KIND=8), INTENT(INOUT) :: Theta(:,:,:,:)
		COMPLEX(KIND=8) :: Gamma(:,:,:), GammaP(:,:,:)
		REAL(KIND=8) :: Lambda(:)
		COMPLEX(KIND=8) :: ThetaTemp(SIZE(Theta,1),localSize,localSize,localSize,SIZE(Gamma,1))
		INTEGER :: alpha, beta, eta, i, j, k, l, d, chi
		chi=SIZE(Gamma,3)
		d=localSize
		DO alpha=1,SIZE(Theta,1)
			DO beta=1,SIZE(Gamma,1)
				DO i=1,d
					DO j=1,d
						DO k=1,d
							ThetaTemp(alpha,i,j,k,beta)=CMPLX(0.0,KIND=8)
							DO eta=1,chi
								ThetaTemp(alpha,i,j,k,beta) = ThetaTemp(alpha,i,j,k,beta) &
															+ Theta(alpha,i,j,eta)*CONJG(Gamma(beta,k,eta)) &
															* Lambda(beta)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
		DO alpha=1,SIZE(GammaP,1)
			DO beta=1,SIZE(Gamma,1)
				DO i=1,d
					DO j=1,d
						Theta(alpha,i,j,beta)=CMPLX(0.0,KIND=8)
						DO k=1,d
							DO eta=1,chi
								Theta(alpha,i,j,beta) = Theta(alpha,i,j,beta) &
													  + Lambda(alpha)*GammaP(alpha,k,eta)*ThetaTemp(eta,i,j,k,beta)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
	END SUBROUTINE ThetaNext

	SUBROUTINE TwoSiteRho(rho2, Theta, Gamma, GammaP, Lambda)
		COMPLEX(KIND=8), INTENT(OUT) :: rho2(:,:)
		COMPLEX(KIND=8), INTENT(IN) :: Theta(:,:,:,:)
		COMPLEX(KIND=8) :: Gamma(:,:,:), GammaP(:,:,:)
		REAL(KIND=8) :: Lambda(:)
		COMPLEX(KIND=8) :: ThetaTemp(SIZE(Theta,1),localSize,localSize,localSize,SIZE(Gamma,1))
		INTEGER :: alpha, beta, eta, i, j, k, l, d
		d=localSize
		DO alpha=1,SIZE(Theta,1)
			DO beta=1,SIZE(Gamma,1)
				DO i=1,d
					DO j=1,d
						DO k=1,d
							ThetaTemp(alpha,i,j,k,beta)=CMPLX(0.0,KIND=8)
							DO eta=1,SIZE(Gamma,3)
								ThetaTemp(alpha,i,j,k,beta) = ThetaTemp(alpha,i,j,k,beta) &
															+ Theta(alpha,i,j,eta)*CONJG(Gamma(beta,k,eta)) &
															* Lambda(beta)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
		DO i=1,d
			DO j=1,d
				DO k=1,d
					DO l=1,d
						rho2((i-1)*d+j,(l-1)*d+k)=CMPLX(0.0,KIND=8)
						DO alpha=1,SIZE(Gamma,1)
							DO beta=1,SIZE(Gamma,3)
								rho2((i-1)*d+j,(l-1)*d+k) = rho2((i-1)*d+j,(l-1)*d+k) &
										+ Lambda(alpha)*GammaP(alpha,i,beta)*ThetaTemp(beta,j,k,l,alpha)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
	END SUBROUTINE TwoSiteRho
	
	FUNCTION TraceMatmul(A,B)
		COMPLEX(KIND=8), INTENT(IN) :: A(:,:), B(:,:)
		COMPLEX(KIND=8) :: TraceMatmul
		INTEGER i, j
		TraceMatmul=CMPLX(0.0,KIND=8)
		DO i=1,SIZE(A,1)
			DO j=1,SIZE(A,2)
				TraceMatmul = TraceMatmul + A(i,j)*B(j,i)
			END DO
		END DO
	END FUNCTION TraceMatmul

	SUBROUTINE TwoSiteExpVal(observable, phaseStat, Op1, Op2, Gammas, Lambdas)
		COMPLEX(KIND=8), INTENT(OUT) :: observable(:,:)
		INTEGER, INTENT(IN) :: phaseStat
		COMPLEX(KIND=8), INTENT(IN) :: Op1(:,:)		
		COMPLEX(KIND=8), INTENT(IN) :: Op2(:,:)
		TYPE(tensor), POINTER :: Gammas(:)
		TYPE(vector), POINTER :: Lambdas(:)
		COMPLEX(KIND=8) :: Theta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3))
		COMPLEX(KIND=8) :: GammaP(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3))
		COMPLEX(KIND=8) :: rho1(localSize,localSize)
		COMPLEX(KIND=8) :: rho2(localSize*localSize,localSize*localSize)
		INTEGER :: l1, l2
		DO l1=1,systemSize
			CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
			observable(l1,l1)=TraceMatmul(Op1,rho1)
		END DO
		DO l2=systemSize,2,(-1)
			CALL ThetaKernal(Theta,Lambdas(l2)%v,Gammas(l2)%t,Lambdas(l2+1)%v)
			DO l1=(l2-1),1,(-1)
				GammaP = Gammas(l1)%t
				! phaseState=2 for fermions
				IF(phaseStat==2) THEN
					CALL OneSiteOp(fermiPhase_op%m,GammaP)
				END IF
				CALL TwoSiteRho(rho2,Theta,Gammas(l1)%t,GammaP,Lambdas(l1)%v)
				observable(l1,l2)=TraceMatmul(Op2,rho2)
				observable(l2,l1)=CONJG(observable(l1,l2))
				CALL ThetaNext(Theta,Gammas(l1)%t,GammaP,Lambdas(l1)%v)
			END DO
		END DO
	END SUBROUTINE TwoSiteExpVal

	SUBROUTINE TotalNumber(number, Gammas, Lambdas)
		REAL(KIND=8), INTENT(OUT) :: number
		TYPE(tensor), POINTER :: Gammas(:)
		TYPE(vector), POINTER :: Lambdas(:)
		TYPE(matrix) :: rho
		INTEGER :: i
		ALLOCATE(rho%m(localSize,localSize))
		number=0.0_8
		DO i=1,systemSize
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
			number=number+REAL(TraceMatmul(n_op%m,rho%m), KIND=8)
		END DO
		DEALLOCATE(rho%m)
	END SUBROUTINE TotalNumber

	SUBROUTINE TotalEnergy(energy, H, Gammas, Lambdas)
		REAL(KIND=8), INTENT(OUT) :: energy
		TYPE(matrix), POINTER :: H(:)
		TYPE(tensor), POINTER :: Gammas(:)
		TYPE(vector), POINTER :: Lambdas(:)
		COMPLEX(KIND=8) :: Theta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3))
		COMPLEX(KIND=8) :: rho2(localSize*localSize,localSize*localSize)
		INTEGER :: l1, l2	
		energy=0.0_8
		DO l2=2,systemSize
			CALL ThetaKernal(Theta,Lambdas(l2)%v,Gammas(l2)%t,Lambdas(l2+1)%v)
			l1=l2-1
			CALL TwoSiteRho(rho2,Theta,Gammas(l1)%t,Gammas(l1)%t,Lambdas(l1)%v)
			energy=energy+REAL(TraceMatmul(H(l1)%m,rho2), KIND=8)
		END DO
	END SUBROUTINE TotalEnergy
	
	SUBROUTINE Qdepletion(depletion, rho, population)
		COMPLEX(KIND=8), INTENT(INOUT) :: rho(:,:)
		REAL(KIND=8), INTENT(IN) :: population
		REAL(KIND=8), INTENT(OUT) :: depletion
		REAL(KIND=8), ALLOCATABLE :: S(:)
		COMPLEX(KIND=8), ALLOCATABLE :: U(:,:), VT(:,:), work(:,:), rwork(:,:)
		INTEGER :: workSize, info
		CHARACTER(1) :: jobu, jobvt
		jobu='A'
		jobvt='A'
		workSize=5*systemSize
		ALLOCATE(S(systemSize))
		ALLOCATE(U(systemSize,systemSize))
		ALLOCATE(VT(systemSize,systemSize))
		ALLOCATE(work(workSize,workSize))
		ALLOCATE(rwork(workSize,workSize))
		CALL ZGESVD(jobu, jobvt, systemSize, systemSize, rho, systemSize, S, U, systemSize, VT, systemSize, &
					work, workSize, rwork, info)
		depletion=1.0_8-S(1)/population
		DEALLOCATE(U,VT,work,rwork,S)
	END SUBROUTINE Qdepletion
	
	SUBROUTINE Egap(gap, rho, densDens, tunnel, population)
		REAL(KIND=8), INTENT(OUT) :: gap
		REAL(KIND=8), INTENT(IN) :: tunnel, population
		COMPLEX(KIND=8), INTENT(INOUT) :: rho(:,:), densDens(:,:)
		REAL(KIND=8) :: k0, pi
		COMPLEX(KIND=8) :: Sk, fsum, Ek, imag
		INTEGER :: i, j
		imag=CMPLX(0.0,1.0,KIND=8)
		pi = 4.0_8*ATAN(1.0_8)
		k0 = pi/(systemSize-1.0_8)
		fsum=CMPLX(0.0,KIND=8)
		DO i=1,(systemSize-1)
		    fsum=fsum + rho(i,i+1) + rho(i+1,i)
		END DO
		Ek=-(tunnel/population)*(COS(k0)-1.0_8)*fsum
		Sk=CMPLX(0.0,KIND=8);
		DO i=1,systemSize
		    DO j=1,systemSize
		        Sk=Sk+EXP(imag*k0*(j-i))*(densDens(i,j)-rho(i,i)*rho(j,j))/population
		    END DO
		END DO
		Sk=Sk+CMPLX(1E-8,KIND=8)
		gap=REAL(Ek/Sk)
	END SUBROUTINE Egap


END MODULE observables_module



