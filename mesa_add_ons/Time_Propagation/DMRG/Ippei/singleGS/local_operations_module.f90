MODULE local_operations_module

	USE system_parameters
	USE MPDtools_module

	IMPLICIT NONE
  
	CHARACTER(1) :: jobu_SVD, jobvt_SVD
	INTEGER :: matrixSizeSM_SVD, workSizeSM_SVD, matrixSizeLG_SVD, workSizeLG_SVD, matrixSize_SVD, workSize_SVD, info_SVD

	REAL(KIND=8), ALLOCATABLE :: rworkSM_SVD(:)
	COMPLEX(KIND=8), ALLOCATABLE :: workSM_SVD(:)
	REAL(KIND=8), ALLOCATABLE :: rworkLG_SVD(:)
	COMPLEX(KIND=8), ALLOCATABLE :: workLG_SVD(:)
	REAL(KIND=8), ALLOCATABLE :: rwork_SVD(:)
	COMPLEX(KIND=8), ALLOCATABLE :: work_SVD(:)


CONTAINS

	SUBROUTINE SVDInit(chi)
		INTEGER, INTENT(IN) :: chi
		jobu_SVD='A'
		jobvt_SVD='A'
		matrixSizeSM_SVD=localSize
		workSizeSM_SVD=5*matrixSizeSM_SVD
		matrixSizeLG_SVD=chi*localSize
		workSizeLG_SVD=5*matrixSizeLG_SVD
		matrixSize_SVD=chi
		workSize_SVD=5*matrixSize_SVD
		ALLOCATE(rworkSM_SVD(workSizeSM_SVD))
		ALLOCATE(workSM_SVD(workSizeSM_SVD))	
		ALLOCATE(rworkLG_SVD(workSizeLG_SVD))
		ALLOCATE(workLG_SVD(workSizeLG_SVD))	
		ALLOCATE(rwork_SVD(workSize_SVD))
		ALLOCATE(work_SVD(workSize_SVD))	
	END SUBROUTINE SVDInit

	SUBROUTINE SVDFinish()	
		DEALLOCATE(rworkSM_SVD)
		DEALLOCATE(workSM_SVD)
		DEALLOCATE(rworkLG_SVD)
		DEALLOCATE(workLG_SVD)
		DEALLOCATE(rwork_SVD)
		DEALLOCATE(work_SVD)
	END SUBROUTINE SVDFinish


	SUBROUTINE OneSiteOp(Op1,Gamma)
		COMPLEX(KIND=8), INTENT(IN) :: Op1(:,:)
		COMPLEX(KIND=8), INTENT(INOUT) :: Gamma(:,:,:)
		INTEGER alpha,chi
		chi = SIZE(Gamma,1);
		DO alpha = 1,chi
			Gamma(alpha,:,:) = MATMUL(Op1,Gamma(alpha,:,:));
		END DO		
	END SUBROUTINE OneSiteOp


	SUBROUTINE FormTheta(Theta, Op2, Gamma1, Lambda1, Gamma2)
		COMPLEX(KIND=8), INTENT(IN) :: Gamma1(:,:,:), Gamma2(:,:,:)
		REAL(KIND=8), INTENT(IN) :: Lambda1(:)
		COMPLEX(KIND=8), INTENT(IN)  :: Op2(:,:)
		COMPLEX(KIND=8), INTENT(OUT) :: Theta(:,:,:,:)
		COMPLEX(KIND=8) :: ThetaTemp(SIZE(Theta,1),SIZE(Theta,2),SIZE(Theta,3),SIZE(Theta,4))
		INTEGER :: chi0,chi1,chi2,alpha,beta,gamma,i,j,k,l
		chi0 = SIZE(Theta,1);
		chi1 = SIZE(Gamma1,3)
		chi2 = SIZE(Theta,4);
		DO alpha=1,chi0
			DO beta=1,chi2
				DO i=1,localSize
					DO j=1,localSize
						ThetaTemp(alpha,i,j,beta) = CMPLX(0.0,KIND=8);					
						DO gamma=1,chi1
							ThetaTemp(alpha,i,j,beta)=ThetaTemp(alpha,i,j,beta) & 
										         + Gamma1(alpha,i,gamma)*Lambda1(gamma)*Gamma2(gamma,j,beta)
						END DO
					END DO
				END DO
			END DO
		END DO
		DO alpha=1,chi0
			DO beta=1,chi2
				DO i=1,localSize
					DO j=1,localSize
						Theta(alpha,i,j,beta) = CMPLX(0.0,KIND=8);					
						DO k=1,localSize
							DO l=1,localSize
								Theta(alpha,i,j,beta)=Theta(alpha,i,j,beta) & 
										             + Op2((i-1)*localSize+j,(k-1)*localSize+l) &
										             * ThetaTemp(alpha,k,l,beta)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
	END SUBROUTINE FormTheta

	SUBROUTINE FormRho(rho,Theta,Lambda0,Lambda2)
		REAL(KIND=8), INTENT(IN) :: Lambda0(:), Lambda2(:)
		COMPLEX(KIND=8), INTENT(IN) :: Theta(:,:,:,:)
		COMPLEX(KIND=8), INTENT(OUT) :: rho(:,:)
		INTEGER :: chi0,chi2,alpha,beta,gamma,i,j,k
		chi0 = SIZE(Theta,1);
		chi2 = SIZE(Theta,4);
		DO alpha=1,chi0
			DO beta=1,chi0
				DO i=1,localSize
					DO j=1,localSize
							rho((i-1)*chi0+alpha,(j-1)*chi0+beta)=CMPLX(0.0,KIND=8)
						DO gamma=1,chi2
							DO k=1,localSize
								rho((i-1)*chi0+alpha,(j-1)*chi0+beta) &
									= rho((i-1)*chi0+alpha,(j-1)*chi0+beta) &
									+ Lambda0(alpha)*Theta(alpha,i,k,gamma)*Lambda2(gamma) &
									* Lambda2(gamma)*CONJG(Theta(beta,j,k,gamma))*Lambda0(beta)
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
	END SUBROUTINE FormRho
	
	SUBROUTINE FormLambda1(Lambda1,S,chi1)
		INTEGER, INTENT(IN) :: chi1
		REAL(KIND=8), INTENT(INOUT) :: Lambda1(:)
		REAL(KIND=8), INTENT(IN) :: S(:)
		INTEGER :: alpha
		DO alpha=1,chi1
			Lambda1(alpha)=SQRT(S(alpha))
			IF(ABS(Lambda1(alpha))<1E-10) THEN
				Lambda1(alpha)=0.0_8
			END IF
		END DO
	END SUBROUTINE FormLambda1

	SUBROUTINE FormGamma1(Lambda0,Gamma1,U,chi0,chi1)
		COMPLEX(KIND=8), INTENT(INOUT) :: Gamma1(:,:,:)
		REAL(KIND=8), INTENT(IN) :: Lambda0(:)
		COMPLEX(KIND=8), INTENT(IN) :: U(:,:)
		INTEGER, INTENT(IN) :: chi0,chi1
		INTEGER :: alpha, beta, i
		DO beta=1,chi1
			DO alpha=1,chi0
				DO i=1,localSize
					IF(ABS(Lambda0(alpha))>1E-10) THEN
						Gamma1(alpha,i,beta)=U((i-1)*chi0+alpha,beta)/Lambda0(alpha)
					ELSE
						Gamma1(alpha,i,beta)=CMPLX(0.0,KIND=8)
					END IF
				END DO
			END DO
		END DO	
	END SUBROUTINE FormGamma1

	SUBROUTINE FormGamma2(Lambda0,Lambda1,Gamma2,U,Theta,chi0,chi2,chi1)
		COMPLEX(KIND=8), INTENT(INOUT) :: Gamma2(:,:,:)
		REAL(KIND=8), INTENT(IN) :: Lambda0(:), Lambda1(:)
		COMPLEX(KIND=8), INTENT(IN) :: U(:,:)
		COMPLEX(KIND=8), INTENT(IN) :: Theta(:,:,:,:)
		INTEGER, INTENT(IN) :: chi0,chi2,chi1
		INTEGER :: alpha, beta, gamma, i, j
		DO alpha=1,chi1
			DO beta=1,chi2
				DO i=1,localSize
					Gamma2(alpha,i,beta)=CMPLX(0.0,KIND=8)
					IF(ABS(Lambda1(alpha))>1E-8) THEN
						DO gamma=1,chi0
							DO j=1,localSize
								Gamma2(alpha,i,beta) = Gamma2(alpha,i,beta) &
			  						+CONJG(U((j-1)*chi0+gamma,alpha))/Lambda1(alpha) &
									*Lambda0(gamma)*Theta(gamma,j,i,beta)
							END DO
						END DO
					END IF
				END DO
			END DO
		END DO
	END SUBROUTINE FormGamma2
	
	SUBROUTINE CanonicalForm(Lambda0,Gamma1,Lambda1,Gamma2,Lambda2)
		COMPLEX(KIND=8), INTENT(INOUT) :: Gamma1(:,:,:), Gamma2(:,:,:)
		REAL(KIND=8), INTENT(INOUT) :: Lambda1(:), Lambda0(:), Lambda2(:)
		COMPLEX(KIND=8) :: M(SIZE(Gamma1,3),SIZE(Gamma1,3)), X(SIZE(Gamma1,3),SIZE(Gamma1,3)), Y(SIZE(Gamma1,3),SIZE(Gamma1,3))
		COMPLEX(KIND=8) :: U(SIZE(Gamma1,3),SIZE(Gamma1,3)), VT(SIZE(Gamma1,3),SIZE(Gamma1,3))
		COMPLEX(KIND=8) :: Gamma1Temp(SIZE(Gamma1,1),SIZE(Gamma1,2),SIZE(Gamma1,3))
		COMPLEX(KIND=8) :: Gamma2Temp(SIZE(Gamma2,1),SIZE(Gamma2,2),SIZE(Gamma2,3))
		REAL(KIND=8) :: S(SIZE(Gamma1,3))
		REAL(KIND=8) :: normsq
		INTEGER :: i,j,k,l,alpha,beta,gamma,chi0,chi1,chi2
		chi0=SIZE(Gamma1,1)
		chi1=SIZE(Gamma1,3)
		chi2=SIZE(Gamma2,3)
		DO alpha=1,chi1
			DO beta=1,chi1
				M(alpha,beta)=CMPLX(0.0,KIND=8)
				DO i=1,localSize
					DO gamma=1,chi0
						M(alpha,beta)=M(alpha,beta)+CONJG(Gamma1(gamma,i,alpha))*Lambda0(gamma)*Lambda0(gamma)*Gamma1(gamma,i,beta)
					END DO
				END DO
			END DO
		END DO
		CALL SVD(M, U, S, VT)
		DO alpha=1,chi0
			DO beta=1,chi1
				DO i=1,localSize
					Gamma1Temp(alpha,i,beta) = CMPLX(0.0,KIND=8)
					IF(S(beta)>1E-8) THEN
						DO gamma=1,chi1
							Gamma1Temp(alpha,i,beta) = Gamma1Temp(alpha,i,beta) &
											     	 + Gamma1(alpha,i,gamma)*U(gamma,beta)*((1.0_8)/SQRT(S(beta)))
						END DO
					END IF
				END DO
			END DO
		END DO
		DO alpha=1,chi1
			DO beta=1,chi1
				X(alpha,beta)=SQRT(S(alpha))*CONJG(U(beta,alpha))
			END DO
		END DO
		DO alpha=1,chi1
			DO beta=1,chi1
				M(alpha,beta)=CMPLX(0.0,KIND=8)
				DO i=1,localSize
					DO gamma=1,chi2
						M(alpha,beta)=M(alpha,beta)+Gamma2(alpha,i,gamma)*Lambda2(gamma)*Lambda2(gamma)*CONJG(Gamma2(beta,i,gamma))
					END DO
				END DO
			END DO
		END DO
		CALL SVD(M, U, S, VT)
		DO alpha=1,chi1
			DO beta=1,chi2
				DO i=1,localSize
					Gamma2Temp(alpha,i,beta) = CMPLX(0.0,KIND=8)
					IF(S(alpha)>1E-8) THEN
						DO gamma=1,chi1
							Gamma2Temp(alpha,i,beta) = Gamma2Temp(alpha,i,beta) &
												 	+ ((1.0_8)/SQRT(S(alpha)))*CONJG(U(gamma,alpha))*Gamma2(gamma,i,beta)
						END DO
					END IF
				END DO
			END DO
		END DO
		DO alpha=1,chi1
			DO beta=1,chi1
				Y(alpha,beta)=Lambda1(alpha)*U(alpha,beta)*SQRT(S(beta))
			END DO
		END DO
		M=MATMUL(X,Y)
		CALL SVD(M, U, S, VT)
		normsq = 0.0_8
		DO alpha=1,chi1
			normsq = normsq + S(alpha)*S(alpha)
		END DO
		Lambda1 = S/DSQRT(normsq)
		DO alpha=1,chi0
			DO beta=1,chi1
				DO i=1,localSize
					Gamma1(alpha,i,beta)=CMPLX(0.0,KIND=8)
					DO gamma=1,chi1
						Gamma1(alpha,i,beta) = Gamma1(alpha,i,beta)+Gamma1Temp(alpha,i,gamma)*U(gamma,beta)
					END DO
				END DO
			END DO
		END DO
		DO alpha=1,chi1
			DO beta=1,chi2
				DO i=1,localSize
					Gamma2(alpha,i,beta) = CMPLX(0.0,KIND=8)
					DO gamma=1,chi1
						Gamma2(alpha,i,beta) = Gamma2(alpha,i,beta) + VT(alpha,gamma)*Gamma2Temp(gamma,i,beta)
					END DO
				END DO
			END DO
		END DO
	END SUBROUTINE CanonicalForm

	SUBROUTINE TwoSiteOp(link,Op2,Gammas,Lambdas)
		INTEGER, INTENT(IN) :: link
		TYPE(tensor), POINTER :: Gammas(:)
		TYPE(vector), POINTER :: Lambdas(:)
		COMPLEX(KIND=8), INTENT(IN) :: Op2(:,:)
		COMPLEX(KIND=8) :: Theta(SIZE(Gammas(link)%t,1),SIZE(Gammas(link)%t,2), &
								 SIZE(Gammas(link+1)%t,2),SIZE(Gammas(link+1)%t,3))
		COMPLEX(KIND=8) :: rho(localSize*SIZE(Gammas(link)%t,1),localSize*SIZE(Gammas(link)%t,1))
		COMPLEX(KIND=8) :: U(localSize*SIZE(Gammas(link)%t,1),localSize*SIZE(Gammas(link)%t,1))
		REAL(KIND=8) :: S(localSize*SIZE(Gammas(link)%t,1))
		INTEGER :: alpha, beta, gamma, i, j, chi1, chi0, chi2, chi
		chi0 = SIZE(Gammas(link)%t,1)
		chi = SIZE(Gammas(1)%t,3)
		chi2 = SIZE(Gammas(link+1)%t,3)
		IF(link/=1) THEN
			chi1=chi
		ELSE
			IF(localSize<=chi) THEN
				chi1=localSize
			ELSE
				chi1=chi
			END IF
		END IF
		CALL FormTheta(Theta,Op2,Gammas(link)%t,Lambdas(link+1)%v,Gammas(link+1)%t)		
		CALL FormRho(rho,Theta,Lambdas(link)%v,Lambdas(link+2)%v)		
		CALL SVDTruncation(link,rho,S,U)
		CALL FormLambda1(Lambdas(link+1)%v,S,chi1)
		CALL FormGamma1(Lambdas(link)%v,Gammas(link)%t,U,chi0,chi1)
		CALL FormGamma2(Lambdas(link)%v,Lambdas(link+1)%v,Gammas(link+1)%t,U,Theta,chi0,chi2,chi1)
!		CALL CanonicalForm(Lambdas(link)%v,Gammas(link)%t,Lambdas(link+1)%v,Gammas(link+1)%t,Lambdas(link+2)%v)
	END SUBROUTINE TwoSiteOp
	
	
	SUBROUTINE SVDTruncation(link,MatrixIn, S, U)
		INTEGER, INTENT(IN) :: link
		COMPLEX(KIND=8), INTENT(INOUT) :: MatrixIn(:,:)       
		REAL(KIND=8), INTENT(OUT) :: S(:)
		COMPLEX(KIND=8), INTENT(OUT) :: U(:,:)
		COMPLEX(KIND=8) :: VT(SIZE(U,1),SIZE(U,2))
		IF(link/=1) THEN
			CALL ZGESVD(jobu_SVD, jobvt_SVD, matrixSizeLG_SVD, matrixSizeLG_SVD, MatrixIn, matrixSizeLG_SVD, S, U, & 
					matrixSizeLG_SVD, VT, matrixSizeLG_SVD, workLG_SVD, workSizeLG_SVD, rworkLG_SVD, info_SVD)
		ELSE
			CALL ZGESVD(jobu_SVD, jobvt_SVD, matrixSizeSM_SVD, matrixSizeSM_SVD, MatrixIn, matrixSizeSM_SVD, S, U, & 
					matrixSizeSM_SVD, VT, matrixSizeSM_SVD, workSM_SVD, workSizeSM_SVD, rworkSM_SVD, info_SVD)
		END IF
	END SUBROUTINE SVDTruncation


	SUBROUTINE SVD(MatrixIn, U, S, VT)
		COMPLEX(KIND=8), INTENT(INOUT) :: MatrixIn(:,:)       
		REAL(KIND=8), INTENT(OUT) :: S(:)
		COMPLEX(KIND=8), INTENT(OUT) :: U(:,:), VT(:,:)
		CALL ZGESVD(jobu_SVD, jobvt_SVD, matrixSize_SVD, matrixSize_SVD, MatrixIn, matrixSize_SVD, S, U, & 
					matrixSize_SVD, VT, matrixSize_SVD, work_SVD, workSize_SVD, rwork_SVD, info_SVD)
	END SUBROUTINE SVD

	



END MODULE local_operations_module













