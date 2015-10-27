!                                           
  MODULE Time_Propagation_Subroutines_Module
  USE Data
  USE Derived_Types
  USE Matrix_Print
  IMPLICIT NONE
!***********************************************************************                                     
!***********************************************************************                                     
!                          Explicit Interfaces                                                               
!***********************************************************************                                    
!                                                                                                           
                            INTERFACE Compute_Energies
                       MODULE PROCEDURE Compute_Energy_Vector_Ham_d,       &
                                        Compute_Energy_Matrix_Ham_d,       &
                                        Compute_Energy_Vector_Ham_z,       &
                                        Compute_Energy_Vector_Exp_d 
                            END INTERFACE Compute_Energies   
!                                                                                                           
                            INTERFACE Initialize_Imaginary_Time
                       MODULE PROCEDURE Initialize_Imaginary_Time_Vector,            &
                                        Initialize_Imaginary_Time_Matrix
                            END INTERFACE Initialize_Imaginary_Time    

!
!*******************************************************************************
                                Contains
!*******************************************************************************
  SUBROUTINE Initialize_Imaginary_Time_Vector (Solution)
  IMPLICIT NONE
  Real(idp),DIMENSION(:)                       :: Solution
  Call Random_Number(Solution)
!*******************************************************************************
  END SUBROUTINE  Initialize_Imaginary_Time_Vector
!*******************************************************************************
  SUBROUTINE Initialize_Imaginary_Time_Matrix (Solution)
  IMPLICIT NONE
  Real(idp),DIMENSION(:,:)                     :: Solution
  Call Random_Number(Solution)
!*******************************************************************************
  END SUBROUTINE  Initialize_Imaginary_Time_Matrix
!*******************************************************************************
  SUBROUTINE Compute_Energy_Vector_Ham_d (hamiltonian, Solution, Mul)
  IMPLICIT NONE
  TYPE(ENE_HAM)                                :: hamiltonian
  REAL(idp),DIMENSION(:)                       :: Solution
  REAL(idp),DIMENSION(:)                       :: MUL
  REAL(idp)                                    :: sdot
  INTEGER                                      :: i
  MUL(1)         = D(1) * Solution(1) + E(1) * Solution(2)
  Do i = 2, M_Size-1
     MUL(i)  = E(i-1) * Solution(i-1) + D(i) * Solution(i) + E(i+1) * Solution(i+1)
  END DO
  i = M_Size
  MUL(i)  = E(i-1) * Solution(i-1) + D(i) * Solution(i) 
  Energy = sdot(M_Size, Solution, 1, MUL, 1)
  write(iout,*) '          The expectation of the Hamiltonian = ', Energy
!*******************************************************************************
  END SUBROUTINE Compute_Energy_Vector_Ham_d
!*******************************************************************************
  SUBROUTINE Compute_Energy_Vector_Ham_z (hamiltonian, Solution, Mul)
  IMPLICIT NONE
  TYPE(ENE_HAM)                                   :: hamiltonian
  COMPLEX(idp),DIMENSION(:)                       :: Solution
  COMPLEX(idp),DIMENSION(:)                       :: MUL
  COMPLEX(idp)                                    :: cdotc
  INTEGER                                         :: i
  MUL(1)         = D(1) * Solution(1) + E(1) * Solution(2)
  Do i = 2, M_Size-1
     MUL(i)  = E(i-1) * Solution(i-1) + D(i) * Solution(i) + E(i+1) * Solution(i+1)
  END DO
  i = M_Size
  MUL(i)  = E(i-1) * Solution(i-1) + D(i) * Solution(i) 
  Energy = cdotc(M_Size, Solution, 1, MUL, 1)
  write(iout,*) '          The expectation of the Hamiltonian = ', Energy
!*******************************************************************************
  END SUBROUTINE Compute_Energy_Vector_Ham_z
!*******************************************************************************
  SUBROUTINE Compute_Energy_Matrix_Ham_d (hamiltonian, Solution, MUL, Ham, Eigval, Eigvec)
  IMPLICIT NONE
  TYPE(ENE_HAM)                                :: hamiltonian
  REAL(idp),DIMENSION(:,:)                     :: Solution
  REAL(idp),DIMENSION(:,:)                     :: Ham
  REAL(idp),DIMENSION(:,:)                     :: Eigvec
  REAL(idp),DIMENSION(:)                       :: Eigval
  REAL(idp),DIMENSION(:)                       :: MUL
  REAL(idp)                                    :: sdot
  REAL(idp)                                    :: overlap
  INTEGER                                      :: i
  INTEGER                                      :: j
!
  DO i = 1, M_Eig
     MUL(1) = D(1) * Solution(1,i) + E(1) * Solution(2,i)
     DO j = 2, M_Size-1
        MUL(j)  = E(j-1) * Solution(j-1,i) + D(j) * Solution(j-1,i) + E(j+1) * Solution(j+1,i)
     END DO
     j = M_Size
     MUL(j)  = E(j-1) * Solution(j-1,i) + D(j) * Solution(j,i)    
     DO j = 1, i
        Ham(j,i) =sdot(M_Size,Solution(1,j),1,MUL,1)
        Ham(i,j) = Ham(j,i)
     END DO
  END DO      
!  Call dsyevx('v','i','u', m_size,Ham, m_size, vl, vu, il, iu, abstol, m_eig, &
!               Eigval, Eigvec, m_size, work_d, lwork, iwork, ifail, info)
!*******************************************************************************      
  END SUBROUTINE Compute_Energy_Matrix_Ham_d
!*******************************************************************************
!*******************************************************************************
  SUBROUTINE Compute_Energy_Vector_Exp_d (exponential, Solution, Previous)
  IMPLICIT NONE
  TYPE(ENE_EXP)                                :: exponential
  REAL(idp),DIMENSION(:)                       :: Solution
  REAL(idp),DIMENSION(:)                       :: Previous
  REAL(idp)                                    :: sdot
  REAL(idp)                                    :: overlap
  INTEGER                                      :: i
  INTEGER                                      :: j
  overlap = sdot(m_size,Solution,1,Solution,1)
  Energy  = sdot(m_size,Solution,1,Previous,1) / overlap
  Energy = log(Energy)/ delta_t
  write(iout,*) '          The expectation of the Hamiltonian = ', Energy
!
!*******************************************************************************      
  END SUBROUTINE Compute_Energy_Vector_Exp_d
!*******************************************************************************
!*******************************************************************************
  SUBROUTINE Comp_E(t)
  IMPLICIT NONE
  INTEGER                                       :: i
  REAL(idp)                                     :: t
  E_t  =  0.0d0
  IF ( t <= Pulse_Time) THEN
       IF (pulse == 'SQUARE' .OR. pulse == 'Square' .OR. pulse == 'square') THEN
           E_t   = E_0
       ELSE IF (pulse == 'SMOOTH' .OR. pulse == 'Smooth' .OR. pulse =='smooth') THEN
           E_t   = E_0 * sin(PI * t / Pulse_Time)**2
       END IF
  ELSE
        write(iout,*) '         Pulse is over.  Time is = ', t
        E_flag = .true.
  END IF
!********************************************************************************
  END SUBROUTINE Comp_E
!********************************************************************************
  SUBROUTINE Comp_A (A, t_l, t_u, pt, wt, WORK)
  IMPLICIT NONE
  INTEGER                                       :: i
  REAL(idp)                                     :: A
  REAL(idp)                                     :: alpha
  REAL(idp)                                     :: beta
  REAL(idp), DIMENSION(2)                       :: endpoints
  REAL(idp)                                     :: t_l
  REAL(idp)                                     :: t_u
  REAL(idp), DIMENSION(Quad_Size)               :: pt
  REAL(idp), DIMENSION(Quad_Size)               :: wt
  REAL(idp), DIMENSION(Quad_Size)               :: WORK
  IF ( t_l <= Pulse_Time .or. t_u <= Pulse_Time) THEN
       CALL gaussq('legendre', Quad_Size, alpha, beta, 0, endpoints, WORK, pt, wt)
       CALL cnvtpt(pt, wt, t_l, t_u, Quad_size)
  END IF
  DO i =1, Quad_Size
     Call Comp_E(pt(i))
     A = A - wt(i) *E_t*sin(Omega*pt(i) + phase)
  END DO
!********************************************************************************
  END SUBROUTINE Comp_A
!********************************************************************************
  SUBROUTINE cnvtpt(pt,wt,t_l,t_u,Quad_size)
  IMPLICIT NONE
  INTEGER                                :: Quad_Size
  REAL(idp), DIMENSION(Quad_Size)        :: pt
  REAL(idp), DIMENSION(Quad_Size)        :: wt
  REAL(idp)                              :: t_l
  REAL(idp)                              :: t_u
  REAL(idp)                              :: f1
  REAL(idp)                              :: f2
  f1 = ( t_u - t_l )*.5D0
  f2 = ( t_l + t_u )*.5D0
  pt =  f1*pt + f2
  wt = wt*f1
!******************************************************************************** 
  END SUBROUTINE cnvtpt
!********************************************************************************
  SUBROUTINE TD_Pot(t)
  IMPLICIT NONE
  REAL(idp)                                  :: t
  REAL(idp)                                  :: fac
!
  fac = E_t *sin(Omega*t + phase)
  V_t(1:m_size) = - fac * x(1:m_size)
!*******************************************************************************            
  END SUBROUTINE TD_Pot
!********************************************************************************
!*******************************************************************************
  SUBROUTINE Get_Probabilities (Solution,prob)
  IMPLICIT NONE
  TYPE(REAL_VECTOR)                               :: type_real_vector
  COMPLEX(idp),DIMENSION(:)                       :: Solution
  COMPLEX(idp)                                    :: temp
  REAL(idp),DIMENSION(:)                          :: prob
  INTEGER                                         :: i
  INTEGER                                         :: j
  DO i = 1, Number_of_Eigenvectors
     temp = (0.d0,0.d0)
     DO j = 1, m_size
        temp = temp+ Solution(j) * Z(j,i)
     END DO
     prob(i) = conjg(temp)*temp
  END DO
  Call Print_Matrix(type_real_vector,prob,title='Probabilities',frmt='fr')
!*******************************************************************************
  END SUBROUTINE Get_Probabilities
!*******************************************************************************
  END MODULE Time_Propagation_Subroutines_Module
!********************************************************************************
