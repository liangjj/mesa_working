MODULE MPDtools_module 

	USE system_parameters_module
	USE Matrix_Exponential_module

	IMPLICIT NONE

! *** DERIVED TYPES ***
  	
	TYPE vector
		REAL(KIND=8), POINTER :: v(:)
	END TYPE vector
	
	TYPE matrix
		COMPLEX(KIND=8), POINTER :: m(:,:)
	END TYPE matrix

	TYPE tensor
		COMPLEX(KIND=8), POINTER :: t(:,:,:)
	END TYPE tensor


! *** Hubbard operators ***
	TYPE(vector) :: x_op
	TYPE(matrix) :: a_op
	TYPE(matrix) :: adag_op
	TYPE(matrix) :: n_op
	TYPE(matrix) :: nsq_op
	TYPE(matrix) :: nxone_op
	TYPE(matrix) :: onexn_op
	TYPE(matrix) :: nsqxone_op
	TYPE(matrix) :: onexnsq_op
	TYPE(matrix) :: adagxa_op
	TYPE(matrix) :: t_op
	TYPE(matrix) :: nxn_op
	TYPE(matrix) :: one_op
	TYPE(matrix) :: onexone_op
	TYPE(matrix) :: fermiPhase_op 


CONTAINS

	SUBROUTINE TensorProduct(AxB,A,B)
		COMPLEX(KIND=8), INTENT(IN) :: A(:,:), B(:,:)
		COMPLEX(KIND=8), INTENT(OUT) :: AxB(:,:)
		INTEGER i,j,k,l, dA, dB
		dA = SIZE(A,1)
		dB = SIZE(B,1)
		DO i=1,dA
			DO j=1,dB
				DO k=1,dA
					DO l=1,dB
						AxB((i-1)*dB+j,(k-1)*dB+l)=A(i,k)*B(j,l)
					END DO
				END DO
			END DO
		END DO
	END SUBROUTINE

	SUBROUTINE CreateFieldOps()
		INTEGER :: i,j,d, dsq
		d=localSize
		dsq=d*d
		ALLOCATE(x_op%v(systemSize))
		ALLOCATE(a_op%m(d,d))
		ALLOCATE(adag_op%m(d,d))
		ALLOCATE(n_op%m(d,d))
		ALLOCATE(nsq_op%m(d,d))
		ALLOCATE(nxone_op%m(dsq,dsq))
		ALLOCATE(onexn_op%m(dsq,dsq))
		ALLOCATE(nsqxone_op%m(dsq,dsq))
		ALLOCATE(onexnsq_op%m(dsq,dsq))
		ALLOCATE(adagxa_op%m(dsq,dsq))
		ALLOCATE(t_op%m(dsq,dsq))
		ALLOCATE(nxn_op%m(dsq,dsq))
		ALLOCATE(one_op%m(d,d))
		ALLOCATE(fermiPhase_op%m(d,d))
		ALLOCATE(onexone_op%m(dsq,dsq))
		a_op%m=CMPLX(0.0,KIND=8)
		adag_op%m=CMPLX(0.0,KIND=8)
		one_op%m=CMPLX(0.0,KIND=8)
		fermiPhase_op%m=CMPLX(0.0,KIND=8)
		DO i=1, systemSize
			x_op%v(i) = i-(systemSize+1.0_8)/2.0_8
		END DO
		DO j=1,(d-1)
			one_op%m(j,j)=CMPLX(1.0,KIND=8)
			fermiPhase_op%m(j,j)=CMPLX((-1.0)**(j-1),KIND=8)
			a_op%m(j,j+1)=CMPLX(SQRT(j*1.0_8),KIND=8)
		END DO
		one_op%m(d,d)=CMPLX(1.0,KIND=8)
		fermiPhase_op%m(d,d)=CMPLX((-1.0)**(d-1),KIND=8)
		adag_op%m = Transpose(a_op%m)
		n_op%m = MATMUL(adag_op%m,a_op%m)
		nsq_op%m = MATMUL(n_op%m,n_op%m)
		CALL TensorProduct(adagxa_op%m,adag_op%m,a_op%m)
		CALL TensorProduct(nxn_op%m,n_op%m,n_op%m)
		CALL TensorProduct(nxone_op%m,n_op%m,one_op%m)
		CALL TensorProduct(onexn_op%m,one_op%m,n_op%m)
		CALL TensorProduct(nsqxone_op%m,nsq_op%m,one_op%m)
		CALL TensorProduct(onexnsq_op%m,one_op%m,nsq_op%m)
		CALL TensorProduct(onexone_op%m,one_op%m,one_op%m)
		t_op%m = adagxa_op%m + Transpose(adagxa_op%m);
	END SUBROUTINE CreateFieldOps

	SUBROUTINE DestroyFieldOps()
		DEALLOCATE(x_op%v)
		DEALLOCATE(a_op%m,adag_op%m,n_op%m,nsq_op%m,one_op%m,fermiPhase_op%m)
		DEALLOCATE(adagxa_op%m, t_op%m, onexone_op%m, nxone_op%m, onexn_op%m, nsqxone_op%m, onexnsq_op%m, nxn_op%m)
	END SUBROUTINE DestroyFieldOps

	SUBROUTINE AllocateGamLam(Gammas, Lambdas, chi)
		TYPE(tensor), POINTER :: Gammas(:)
		TYPE(vector), POINTER :: Lambdas(:)
		INTEGER, INTENT(IN) :: chi
		INTEGER :: i		
		ALLOCATE(Gammas(systemSize))
		ALLOCATE(Lambdas(systemSize+1))
		ALLOCATE(Lambdas(1)%v(1))
		Lambdas(1)%v=0.0_8
		ALLOCATE(Gammas(1)%t(1,localSize,chi));
		Gammas(1)%t=CMPLX(0.0,KIND=8)
		DO	i=2,(systemSize-1)
			ALLOCATE(Gammas(i)%t(chi,localSize,chi))
			Gammas(i)%t=CMPLX(0.0,KIND=8)
			ALLOCATE(Lambdas(i)%v(chi))
			Lambdas(i)%v=CMPLX(0.0,KIND=8)
		END DO
		ALLOCATE(Lambdas(systemSize)%v(chi))
		Lambdas(systemSize)%v=0.0_8
		ALLOCATE(Lambdas(systemSize+1)%v(1))
		Lambdas(systemSize+1)%v=0.0_8
		ALLOCATE(Gammas(systemSize)%t(chi,localSize,1))
		Gammas(systemSize)%t=CMPLX(0.0,KIND=8)
	END SUBROUTINE AllocateGamLam

	SUBROUTINE CopyGamLam(GammasCopy, LambdasCopy, GammasOrig, LambdasOrig)
		TYPE(tensor), POINTER :: GammasCopy(:), GammasOrig(:)
		TYPE(vector), POINTER :: LambdasCopy(:), LambdasOrig(:)
		INTEGER :: i, alpha, j, beta, chimin
		chimin=MIN(SIZE(GammasCopy(1)%t,3),SIZE(GammasOrig(1)%t,3))
		LambdasCopy(1)%v(1)=1.0_8
		LambdasCopy(systemSize+1)%v(1)=1.0_8
		DO i=2,systemSize
			DO alpha=1,chimin
				LambdasCopy(i)%v(alpha)=LambdasOrig(i)%v(alpha)
			END DO
		END DO
		DO alpha=1,chimin
			DO j=1,localSize
				GammasCopy(1)%t(1,j,alpha)=GammasOrig(1)%t(1,j,alpha)
				GammasCopy(systemSize)%t(alpha,j,1)=GammasOrig(systemSize)%t(alpha,j,1)
			END DO
		END DO
		DO i=2,(systemSize-1)
			DO alpha=1,chimin
				DO beta=1,chimin
					DO j=1,localSize
						GammasCopy(i)%t(alpha,j,beta)=GammasOrig(i)%t(alpha,j,beta)
					END DO
				END DO
			END DO
		END DO
	END SUBROUTINE CopyGamLam

	SUBROUTINE DeallocateGamLam(Gammas, Lambdas)
		TYPE(tensor), POINTER :: Gammas(:)
		TYPE(vector), POINTER :: Lambdas(:)	
		INTEGER :: i
		DO	i=1,(systemSize)
			DEALLOCATE(Gammas(i)%t)
			DEALLOCATE(Lambdas(i)%v)
		END DO
		DEALLOCATE(Lambdas(systemSize+1)%v)
		DEALLOCATE(Gammas)
		DEALLOCATE(Lambdas)
	END SUBROUTINE DeallocateGamLam
	
	SUBROUTINE FockState(args, Gammas, Lambdas)
		INTEGER, INTENT(IN) :: args(:)
		TYPE(tensor), POINTER :: Gammas(:)
		TYPE(vector), POINTER :: Lambdas(:)	
		INTEGER :: i, j
		COMPLEX(KIND=8) :: stateFunction
		DO i=1,systemSize
			Gammas(i)%t=CMPLX(0.0, KIND=8)
			Lambdas(i)%v=0.0_8
			Lambdas(i)%v(1)=1.0_8
			DO j=1,localSize
				IF(j==(args(i)+1)) THEN
					Gammas(i)%t(1,j,1)=CMPLX(1.0,KIND=8)
				ELSE
					Gammas(i)%t(1,j,1)=CMPLX(0.0,KIND=8)
				END IF			
			END DO
		END DO
		Lambdas(systemSize+1)%v=0.0_8
		Lambdas(systemSize+1)%v(1)=1.0_8
	END SUBROUTINE FockState

	SUBROUTINE AllStates(Gammas, Lambdas)
		TYPE(tensor), POINTER :: Gammas(:)
		TYPE(vector), POINTER :: Lambdas(:)	
		INTEGER :: i, j
		DO i=1,systemSize
			Gammas(i)%t=CMPLX(0.0, KIND=8)
			Lambdas(i)%v=0.0_8
			Lambdas(i)%v(1)=1.0_8
			DO j=1,localSize
				Gammas(i)%t(1,j,1)=CMPLX((1.0_8)/SQRT(localSize*1.0_8),KIND=8)
			END DO
		END DO
		Lambdas(systemSize+1)%v=0.0_8
		Lambdas(systemSize+1)%v(1)=1.0_8
	END SUBROUTINE AllStates
	
	SUBROUTINE AllocateOps(Ops,numops,opsize)
		TYPE(matrix), POINTER :: Ops(:)
		INTEGER, INTENT(IN) :: numops,opsize
		INTEGER :: i
		ALLOCATE(Ops(numops))
		DO i=1,numops
			ALLOCATE(Ops(i)%m(opsize,opsize))
		END DO
	END SUBROUTINE AllocateOps

	SUBROUTINE DeallocateOps(Ops,numops)
		TYPE(matrix), POINTER :: Ops(:)
		INTEGER, INTENT(IN) :: numops
		INTEGER :: i
		DO i=1,numops
			DEALLOCATE(Ops(i)%m)
		END DO
		DEALLOCATE(Ops)
	END SUBROUTINE DeallocateOps

	SUBROUTINE HamiltonianBoseHubbard(H,trap,tunnel,interaction,mu,x0,lambda)
		TYPE(matrix), POINTER :: H(:)
		REAL(KIND=8), INTENT(IN) :: trap, tunnel, interaction, mu, x0, lambda
		INTEGER :: i, n
		REAL(KIND=8) :: xi, xj, lamtun
		n=systemSize
		lamtun=lambda   !*tunnel
		DO i=1,(n-1)
			xi=x_op%v(i) + x0;
			xj=x_op%v(i+1) + x0;
			H(i)%m = ((trap*xi*xi + 2.0_8*lamtun*COS(3.14159265358979*i) - mu - interaction/2.0_8)*nxone_op%m &
				   + interaction/2.0_8*nsqxone_op%m)/2.0_8 &
				   + ((trap*xj*xj + 2.0_8*lamtun*COS(3.14159265358979*(i+1)) - mu - interaction/2.0_8)*onexn_op%m &
				   + interaction/2.0_8*onexnsq_op%m)/2.0_8 &
				   - tunnel*t_op%m
		END DO
		xi=x_op%v(1) + x0
		H(1)%m = H(1)%m + ((trap*xi*xi + 2.0_8*lamtun*COS(3.14159265358979) - mu - interaction/2.0_8)*nxone_op%m &
		       + interaction/2.0_8*nsqxone_op%m)/2.0_8
		xi=x_op%v(n) + x0
		H(n-1)%m = H(n-1)%m + ((trap*xi*xi + 2.0_8*lamtun*COS(3.14159265358979*n) - mu - interaction/2.0_8)*onexn_op%m &
		         + interaction/2.0_8*onexnsq_op%m)/2.0_8
	END SUBROUTINE HamiltonianBoseHubbard

	SUBROUTINE ConstructPropagators(H, U, dtodd, dteven)
		TYPE(matrix), POINTER :: H(:)
		TYPE(matrix), POINTER :: U(:)
		COMPLEX(KIND=8), INTENT(IN) :: dtodd, dteven
		INTEGER :: i
		DO i=1,(systemSize-1),2
			CALL Matrix_Exponential(H(i)%m, U(i)%m, dtodd, localSize*localSize)
		END DO
		DO i=2,(systemSize-1),2
			CALL Matrix_Exponential(H(i)%m, U(i)%m, dteven, localSize*localSize)
		END DO
	END SUBROUTINE ConstructPropagators


END MODULE MPDtools_module














