!
  module SO_Module
  USE Data
  USE Matrix_Print
  USE Time_Propagation_Subroutines_Module
  IMPLICIT NONE
!
!                    Recall That the default matrix operations are;
!                Exp(odd*t/2) * Exp(even*t) * Exp(even*t/2)
!                So the times must be chosen correctly in the propagator
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
!
                            INTERFACE SO_Propagation_Driver                      
                       MODULE PROCEDURE SO_Driver_Length,                  &
                                        SO_Driver_Velocity,                &
                                        SO_Driver_Imag  
                            END INTERFACE SO_Propagation_Driver
!
                            INTERFACE Prop                       
                       MODULE PROCEDURE Prop_len,                          &
                                        Prop_vel,                          &  
                                        Prop_imag
                            END INTERFACE Prop
!
                            INTERFACE Sector                       
                       MODULE PROCEDURE Sector_len,                        &
                                        Sector_vel,                        &  
                                        Sector_imag
                            END INTERFACE Sector
!
                            INTERFACE V_Scale                       
                       MODULE PROCEDURE V_Scale_d,                         &
                                        V_Scale_z,                         &
                                        V_Scale_d_mat                         
                            END INTERFACE V_Scale
!
                            INTERFACE Exp_Off_Diag                       
                       MODULE PROCEDURE Exp_Off_Diag_z,                    &
                                        Exp_Off_Diag_d,                    &
                                        Exp_Off_Diag_d_mat
                            END INTERFACE Exp_Off_Diag                       
!
                            INTERFACE Mat_2_V_2                       
                       MODULE PROCEDURE Mat_2_V_2_z,                       &
                                        Mat_2_V_2_d,                       &
                                        Mat_2_Mat_2_d                       
                            END INTERFACE Mat_2_V_2
!
                            INTERFACE Orthonormalize                       
                       MODULE PROCEDURE Orthonormalize_z,                  &
                                        Orthonormalize_d,                  &
                                        Normalize_d,                       &
                                        Normalize_z
                            END INTERFACE Orthonormalize
!
!____________________________________________________________________________________________!
!
  INTEGER                                                :: trips
  INTEGER                                                :: n_trips
  INTEGER                                                :: mat_size
  INTEGER, DIMENSION(3)                                  :: locate
  INTEGER, DIMENSION(3)                                  :: begin
  REAL(idp), DIMENSION(3)                                :: t_b
  REAL(idp)                                              :: diff
  REAL(idp)                                              :: E_old
!********************************************************************************
                                    Contains
!********************************************************************************
!********************************************************************************
  SUBROUTINE SO_Driver 
  IMPLICIT NONE
  TYPE(SPLIT_OPERATOR)                                :: so
  TYPE(REAL_TIME_SO_LEN)                              :: so_len
  TYPE(REAL_TIME_SO_VEL)                              :: so_vel
  TYPE(IMAG_TIME_SO)                                  :: so_imag
  TYPE(MATRICES), DIMENSION(:),   ALLOCATABLE         :: mat 
  TYPE(COMPLEX_H)                                     :: so_z
  TYPE(REAL_H)                                        :: so_d
  TYPE(COMPLEX_PROP)                                  :: pr_z
  TYPE(REAL_PROP)                                     :: pr_d
  INTEGER                                             :: len
  IF (Method == 'real_time_so_length') THEN
      Call SO_Propagation_Driver(so, so_len, mat, so_z, so_d, pr_z)                       
  ELSE IF (Method == 'real_time_so_velocity') THEN
      Call SO_Propagation_Driver(so, so_vel, mat, so_z, so_d, pr_z)                           
  ELSE IF (Method == 'imaginary_time_so') THEN    
      Call SO_Propagation_Driver(so, so_imag, mat, so_d, pr_d)                       
  ELSE
      Call lnkerr('bad call to SO routines')
  END IF
!********************************************************************************
  END SUBROUTINE SO_Driver 
!********************************************************************************
  SUBROUTINE SO_Driver_length(so, so_len, mat, so_z, so_d, pr_z)
  IMPLICIT NONE
  TYPE(SPLIT_OPERATOR)                           :: so
  TYPE(REAL_TIME_SO_LEN)                         :: so_len
  TYPE(MATRICES), DIMENSION(:), ALLOCATABLE      :: mat 
  TYPE(COMPLEX_H)                                :: so_z
  TYPE(REAL_H)                                   :: so_d
  TYPE(COMPLEX_PROP)                             :: pr_z      
  TYPE(COMPLEX_VECTOR)                           :: type_complex_vector      
  TYPE(ENE_HAM)                                  :: hamiltonian
  INTEGER                                        :: i
!
!
  write(iout,1)
  write(iout,2)
  write(iout,1)
  lwork = 3*m_size - 2
  ALLOCATE( so_d%ham_d(1:m_size), so_d%ham_l(1:m_size), so_d%ham_u(1:m_size), so_z%sol(1:m_size), &
            so_d%eigvec(1:2,1:2), so_d%eigval(1:2), so_z%scr(1:m_size), work_z(1:lwork), rwork(1:lwork) )
!
  ALLOCATE(so%mat(1:M_Size-1))
  DO i = 1, M_Size-1
     ALLOCATE(so%mat(i)%pr_z%prop(1:2,1:2))
  END DO
!
  so_z%sol(1:m_size) = Z(1:m_size,1)
  so_d%ham_d(1:m_size) = D(:)
  so_d%ham_l(:) = E(:)
  so_d%ham_u(:) = E(:)
  call cpu_time (t_i)
  pulse_time  =  delta_t * m_size * 0.5d0
  time = 0.0d0
!
! Set sector times
!
  t_b(1) = delta_t * .5d0  !  note the half time for odd blocks and the potential
  t_b(2) = delta_t 
  begin(1) = 1
  begin(2) = 2
  begin(3) = 1
  locate(1) = 1
  locate(2) = 2
  locate(3) = 1
  n_trips = 3
  mat_size = m_size - 1
  IF ( locate(2) == 1 ) THEN
       n_trips =1
  END IF
  DO i = 1, NO_Time_Steps
     t_l = time + .5d0*delta_t
     IF ( E_flag /= .true. ) THEN
          Call Comp_E (t_l)
          Call TD_Pot(t_l)
     END IF
     Call Prop( so, mat, so_d, pr_z)
     Call V_Scale(so_z%sol, v)
     Call Exp_Off_Diag(so, mat, so_z, pr_z, so_z%sol)
     Call V_Scale(so_z%sol, v)
     time  = time + delta_t
  END DO
  call cpu_time (t_f)
  write (iout,*) '***** Length Solution *****'
  IF ( t_f - t_i < 60 ) THEN
     write (iout,*) 'Time elapsed: ', (t_f - t_i),'seconds'
  ELSE
     write (iout,*) 'Time elapsed: ', (t_f - t_i)/real(60),'minutes'
  END IF
  call Compute_Energies (hamiltonian,so_z%sol,so_z%scr)
  write (iout,*)
  write (iout,*) '"TEST 1"'
  write (iout,*) 'Energy AT E = ', E_0, ': ', Energy
!-------------------------------------------------------------------------!
  IF (Prnt(3) == .true. ) THEN
      Call Print_Matrix(type_complex_vector,so_z%sol,title='SO Length Gauge Solution Vector:', &
                        frmt='fr')
  END IF
  DO i = 1, M_Size - 1
     DEALLOCATE(so%mat(i)%pr_z%prop)
  END DO
  DEALLOCATE(so%mat)
  DEALLOCATE( so_d%ham_d, so_d%ham_l, so_d%ham_u, so_z%sol, so_d%eigvec, so_d%eigval,       &
              so_z%scr, work_z, rwork )
1 Format(//,'********************************************************************************')
2 Format(/,20x,' Beginning Real Time SO Propagation Using Length Gauge')
!*******************************************************************************
  END SUBROUTINE SO_Driver_Length
!********************************************************************************
  SUBROUTINE SO_Driver_Velocity(so, so_vel, mat, so_z, so_d, pr_z)
  IMPLICIT NONE
  TYPE(SPLIT_OPERATOR)                           :: so
  TYPE(REAL_TIME_SO_VEL)                         :: so_vel
  TYPE(MATRICES), DIMENSION(:), ALLOCATABLE      :: mat 
  TYPE(COMPLEX_H)                                :: so_z
  TYPE(REAL_H)                                   :: so_d
  TYPE(COMPLEX_PROP)                             :: pr_z      
  TYPE(COMPLEX_VECTOR)                           :: type_complex_vector      
  TYPE(ENE_HAM)                                  :: hamiltonian
  INTEGER                                        :: i  
!
!
  write(iout,1)
  write(iout,2)
  write(iout,1)
  pulse_time  =  delta_t * m_size * 0.5d0
  quad_size = 4
  lwork = 3*m_size - 2
  ALLOCATE( so_z%ham_d(1:m_size), so_z%ham_l(1:m_size), so_z%ham_u(1:m_size), so_z%sol(1:m_size), &
            so_z%eigvec(1:2,1:2), so_d%eigval(1:2), so_z%scr(1:m_size), work_z(1:lwork), rwork(1:lwork) )
!
  ALLOCATE(so%mat(1:M_Size-1))
  DO i = 1, M_Size-1
     ALLOCATE(so%mat(i)%pr_z%prop(1:2,1:2))
  END DO
  so_z%sol(1:m_size) = Z(1:m_size,1)
  A = 0.d0
  CAll Comp_A(A,0.d0,delta_t*.5d0,pt,wt,rwork)
  call cpu_time (t_i)
  time = 0.0d0
!
! Set sector times
!
  t_b(1) = delta_t * .5d0  !  note the half time for odd blocks and the potential
  t_b(2) = delta_t 
  begin(1) = 1
  begin(2) = 2
  begin(3) = 1
  locate(1) = 1
  locate(2) = 2
  locate(3) = 1
  n_trips = 3
  mat_size = m_size - 1
  IF ( locate(2) == 1 ) THEN
       n_trips =1
  END IF
  DO i =1, No_time_steps
     write(iout,*) '   Time Step = ',i
     A_old = A
     so_z%ham_d(1:M_Size)   = D(1:M_Size) + A*A*0.5d0
     so_z%ham_l(1: M_Size)  = E(1:M_Size) - eye * A * 0.5d0 / Step_Size
     so_z%ham_u(1:M_Size)   = E(1:M_Size) + eye * A * 0.5d0 / Step_Size
     Call Prop(so, mat, so_z, pr_z)
     Call Exp_Off_Diag(so, mat, so_z, pr_z, so_z%sol)
     CALL Comp_A (A, time + delta_t*0.5d0, time + delta_t*1.5d0, pt, wt, rwork)
     time  = time + delta_t
  END DO
  call cpu_time (t_f)
  write (iout,*) '***** Velocity Solution *****'
  IF ( t_f - t_i < 60 ) THEN
     write (iout,*) 'Time elapsed: ', (t_f - t_i),'seconds'
  ELSE
     write (iout,*) 'Time elapsed: ', (t_f - t_i)/real(60),'minutes'
  END IF
  Call Compute_Energies (hamiltonian,so_z%sol,so_z%scr)
  write (iout,*)
  write (iout,*) '"TEST 1"'
  write (iout,*) 'Energy AT E = ', E_0, ': ', Energy
!-------------------------------------------------------------------------!
  IF (Prnt(3) == .true. ) THEN
      Call Print_Matrix(type_complex_vector,so_z%sol,title='SO Velocity Gauge Solution Vector:', &
                        frmt='fr')
  END IF
  DO i = 1, M_Size - 1
     DEALLOCATE(so%mat(i)%pr_z%prop)
  END DO
  DEALLOCATE(so%mat)
  DEALLOCATE( so_z%ham_d, so_z%ham_l, so_z%ham_u, so_z%sol, so_z%eigvec, so_d%eigval,        &
              work_z, rwork )
1 Format(//,'********************************************************************************')
2 Format(/,20x,' Beginning Real Time SO Velocity Gauge Propagation')
!*******************************************************************************
  END SUBROUTINE SO_Driver_Velocity
!********************************************************************************
  SUBROUTINE SO_Driver_Imag(so, so_imag, mat, so_d, pr_d)
  IMPLICIT NONE
  TYPE(SPLIT_OPERATOR)                           :: so
  TYPE(IMAG_TIME_SO)                             :: so_imag
  TYPE(MATRICES), DIMENSION(:), ALLOCATABLE      :: mat 
  TYPE(REAL_H)                                   :: so_d
  TYPE(REAL_PROP)                                :: pr_d      
  TYPE(REAL_VECTOR)                              :: type_real_vector      
  TYPE(ENE_HAM)                                  :: hamiltonian
  TYPE(ENE_EXP)                                  :: exponential
  INTEGER                                        :: i  
  REAL(idp)                                      :: sdot
  REAL(idp)                                      :: overlap
!
!
  write(iout,1)
  write(iout,2)
  write(iout,1)
  lwork = 10*m_size
!
!  Break up the matrix into M_size-1 overlapping (2,2) blocks.
!
  ALLOCATE(so%mat(1:M_Size-1))
  DO i = 1, M_Size-1
     ALLOCATE(so%mat(i)%pr_d%prop(1:2,1:2))
  END DO

  call cpu_time (t_i)
  time = 0.0d0
!
! Set sector times
!
  t_b(1) = delta_t * .5d0  !  note the half time for odd blocks and the potential
  t_b(2) = delta_t 
  begin(1) = 1
  begin(2) = 2
  begin(3) = 1
  locate(1) = 1
  locate(2) = 2
  locate(3) = 1
  n_trips = 3
  mat_size = m_size - 1
  IF ( locate(2) == 1 ) THEN
       n_trips =1
  END IF 
!
  IF (M_Eig == 1 ) THEN
      ALLOCATE( so_d%ham_d(1:m_size), so_d%ham_l(1:m_size), so_d%ham_u(1:m_size), so_d%sol(1:m_size), &
                so_d%sol_v(1:m_size), so_d%eigvec(1:2,1:2), so_d%eigval(1:2), so_d%scr(1:m_size),     &
                work_d(1:lwork) )
      so_d%ham_d(1:m_size) = D(:)
      so_d%ham_l(:) = E(:)
      so_d%ham_u(:) = E(:)
!      Call Initialize_Imaginary_Time(so_d%sol)
!      Z(:,1)=0.0
!      Z(1,1) = 1.d0
      so_d%sol(:) = Z(:,1)
!      so_d%sol = Z(:,1)
!      so_d%sol_v(:) = so_d%sol(:)
!      Call Compute_Energies(exponential,so_d%sol,so_d%sol_v)
      Call Compute_Energies(hamiltonian,so_d%sol,so_d%scr)
      E_old=Energy
      diff = 1.d0
      i  = 0
      DO while ( diff > E_con .and. i <= NO_Time_Steps)
         i = i + 1
         write(iout,*) '   Time Step = ',i
         Call Prop(so, mat, so_d, pr_d)
         Call V_Scale(so_d%sol, v)
         Call Exp_Off_Diag(so, mat, so_d, pr_d, so_d%sol)
         Call V_Scale(so_d%sol, v)
         Call Orthonormalize( so_d%sol, m_size)
!         Call Compute_Energies(exponential, so_d%sol, so_d%sol_v)
         Call Compute_Energies(hamiltonian, so_d%sol, so_d%scr)
         diff = abs(Energy - E_old)
         E_old = Energy
         time  = time + delta_t
!         so_d%sol_v = so_d%sol
      END DO
  ELSE
!      ALLOCATE( so_d%ham_d(1:m_size), so_d%ham_l(1:m_size), so_d%ham_u(1:m_size),                     &
!                so_d%sol_m(1:m_size,1:m_eig), so_d%h_mat(1:m_size,1:m_size), so_d%h_eigval(1:m_eig),  &
!                so_d%h_eigvec(1:m_size,1:m_eig), so_d%eigvec(1:2,1:2), so_d%eigval(1:2),              &
!                so_d%scr(1:m_size), work_d(1:lwork), iwork(1:5*m_size) )
      so_d%ham_d(1:m_size) = D(:)
      so_d%ham_l(:) = E(:)
      so_d%ham_u(:) = E(:)
      Call Initialize_Imaginary_Time(so_d%sol_m)
      Call Compute_Energies(hamiltonian,so_d%sol_m,so_d%scr,so_d%H_mat,so_d%eigval,so_d%eigvec)
  END IF
!

  Write(iout,1)
  Write (iout,*) '         The Converged Energy = ', Energy
  Write(iout,1)
  call cpu_time (t_f)
  IF ( t_f - t_i < 60 ) THEN
     write (iout,*) 'Time elapsed: ', (t_f - t_i),'seconds'
  ELSE
     write (iout,*) 'Time elapsed: ', (t_f - t_i)/real(60),'minutes'
  END IF
!-------------------------------------------------------------------------!
  IF (Prnt(3) == .true. ) THEN
      Call Print_Matrix(type_real_vector,so_d%sol,title='SO Imaginary Time Solution Vector:', &
                        frmt='fr')
  END IF
  DO i = 1, M_Size - 1
     DEALLOCATE(so%mat(i)%pr_d%prop)
  END DO
  DEALLOCATE(so%mat)
  DEALLOCATE( so_d%ham_d, so_d%ham_l, so_d%ham_u, so_d%sol, so_d%eigvec, so_d%eigval,        &
              so_d%scr, work_d )
1 Format(/,'********************************************************************************')
2 Format(/,20x,' Beginning SO Imaginary Time Propagation')
!********************************************************************************
  END SUBROUTINE SO_Driver_Imag
!********************************************************************************
  SUBROUTINE Prop_Vel(so, mat, so_z, pr_z)
  IMPLICIT NONE
  TYPE(SPLIT_OPERATOR)                           :: so
  TYPE(REAL_TIME_SO_VEL)                         :: so_vel
  TYPE(MATRICES), DIMENSION(:), ALLOCATABLE      :: mat
  TYPE(COMPLEX_H)                                :: so_z
  TYPE(COMPLEX_PROP)                             :: pr_z
  INTEGER                                        :: i
  CHARACTER (LEN=5)                              :: itoc 
!
! The matrix size is m_size which means there are m_size - 1 blocks
!
  so_z%ham_d(2:m_size-1) =  so_z%ham_d(2:m_size-1) * .5d0
!
  DO trips = 1, 2
     DO i = begin(trips), M_Size- 1, 2
        so_z%eigvec(1,1) = so_z%ham_d(i)
        so_z%eigvec(2,2) = so_z%ham_d(i+1)
        so_z%eigvec(1,2) = so_z%ham_u(i)
        so_z%eigvec(2,1) = so_z%ham_l(i)
        Call zheev('v','l',2,so_z%eigvec,2,so_z%eigval,work_z,lwork,rwork,info)    
        Call Sector(so%mat(i)%pr_z%prop, so_z%eigval, so_z%eigvec, t_b(trips))
     END DO
  END DO
!
!*******************************************************************************
  END SUBROUTINE  Prop_Vel
!*******************************************************************************
  SUBROUTINE Prop_Len(so, mat, so_d, pr_z)
  IMPLICIT NONE
  TYPE(SPLIT_OPERATOR)                           :: so
  TYPE(MATRICES), DIMENSION(:), ALLOCATABLE      :: mat 
  TYPE(COMPLEX_PROP)                             :: pr_z
  TYPE(REAL_H)                                   :: so_d      
  INTEGER                                        :: i
  CHARACTER (LEN=5)                              :: itoc 
!
!        The matrix size is m_size which means there are m_size - 1 blocks
!
!        Calculate eigenvalues, vectors and propagators of the odd blocks
!
  so_d%ham_d(2:m_size-1) =  so_d%ham_d(2:m_size-1) * .5d0
  DO trips = 1, 2
     DO i = begin(trips), M_Size- 1, 2
        so_d%eigvec(1,1) = so_d%ham_d(i)
        so_d%eigvec(2,2) = so_d%ham_d(i+1)
        so_d%eigvec(1,2) = so_d%ham_u(i)
        so_d%eigvec(2,1) = so_d%ham_l(i)
        Call dsyev('v','l',2,so_d%eigvec,2,so_d%eigval,work_d,lwork,info)    
        Call Sector(so%mat(i)%pr_z%prop, so_d%eigval, so_d%eigvec, t_b(trips))
     END DO
  END DO
!
!********************************************************************************
  END SUBROUTINE Prop_Len
!********************************************************************************
  SUBROUTINE Prop_Imag(so, mat, so_d, pr_d)
  IMPLICIT NONE
  TYPE(SPLIT_OPERATOR)                           :: so
  TYPE(IMAG_TIME_SO)                             :: so_imag
  TYPE(MATRICES), DIMENSION(:), ALLOCATABLE      :: mat       
  TYPE(REAL_PROP)                                :: pr_d
  TYPE(REAL_H)                                   :: so_d
  INTEGER                                        :: i
  CHARACTER (LEN=5)                              :: itoc 
  TYPE(REAL_VECTOR)                              :: type_real_vector      
  TYPE(REAL_MATRIX)                              :: type_real_matrix      
!
  so_d%ham_d(1) = D(1)
  so_d%ham_d(m_size) = D(m_size)
  so_d%ham_d(2:m_size-1) =  D(2:m_size-1) * .5d0
  DO trips = 1, 2
     DO i = begin(trips), M_Size-1, 2
        so_d%eigvec(1,1) = so_d%ham_d(i)
        so_d%eigvec(2,2) = so_d%ham_d(i+1)
        so_d%eigvec(1,2) = so_d%ham_u(i)
        so_d%eigvec(2,1) = so_d%ham_l(i)
        IF(Prnt(4) == .true.) THEN
           write(iout,*) '**********************************************************************'
           write(iout,*) '                        Sector = ',i
           write(iout,*) '**********************************************************************'
           write(iout,*)
           Call Print_Matrix(type_real_matrix,so_d%eigvec,2,2,title='sector Hamiltonian')
       END IF
        Call dsyev('v','u',2,so_d%eigvec,2,so_d%eigval,work_d,lwork,info)    
        IF(Prnt(4) == .true.) THEN
           write(iout,*)
           Call Print_Matrix(type_real_vector,so_d%eigval,frmt='fr',title='sector eigenvalues')
           write(iout,*)
           Call Print_Matrix(type_real_matrix,so_d%eigvec,2,2,title='sector eigenvectors')
        END IF
        Call Sector(so%mat(i)%pr_d%prop, so_d%eigval, so_d%eigvec, t_b(trips))
        IF(Prnt(5) == .true.) THEN
           write(iout,*)
           Call Print_Matrix(type_real_matrix,so%mat(i)%pr_d%prop,2,2,title='sector propagator')
        END IF
     END DO
  END DO
!********************************************************************************
  END SUBROUTINE Prop_Imag
!********************************************************************************
  SUBROUTINE Sector_vel( prop, eigval, eigvec, t)
  IMPLICIT NONE
  COMPLEX(idp), DIMENSION(:,:)                   :: prop
  REAL(idp), DIMENSION(:)                        :: eigval
  COMPLEX(idp), DIMENSION(:,:)                   :: eigvec
  REAL(idp)                                      :: t
  COMPLEX(idp), DIMENSION(2,2)                   :: tmp
  INTEGER                                        :: i
  INTEGER                                        :: j
  INTEGER                                        :: lambda
  prop(:,:) = (0.d0,0.d0)
  DO i = 1, 2
     DO lambda = 1, 2  
        tmp(lambda,i) = eigvec(i,lambda) * exp(-eye*eigval(lambda)*t)
     END DO
  END DO
  DO i = 1, 2
     DO j = 1, 2
        DO lambda = 1 ,2
           prop(i,j) = prop(i,j) + eigvec(i,lambda) * tmp(lambda,j) 
        END DO
     END DO
  END DO
!********************************************************************************
  END SUBROUTINE Sector_vel
!********************************************************************************
!*******************************************************************************
  SUBROUTINE Sector_len( prop, eigval, eigvec, t)
  IMPLICIT NONE
  COMPLEX(idp), DIMENSION(:,:)                   :: prop
  REAL(idp), DIMENSION(:)                        :: eigval
  REAL(idp), DIMENSION(:,:)                      :: eigvec
  REAL(idp)                                      :: t
  REAL(idp), DIMENSION(2,2)                      :: tmp
  INTEGER                                        :: i
  INTEGER                                        :: j
  INTEGER                                        :: lambda
  prop(:,:) = (0.d0,0.d0)
  DO i = 1, 2
     DO lambda = 1, 2  
        tmp(lambda,i) = eigvec(i,lambda) * exp(-eye*eigval(lambda)*t)
     END DO
  END DO
  DO i = 1, 2
     DO j = 1, 2
        DO lambda = 1 ,2
           prop(i,j) = prop(i,j) + eigvec(i,lambda) * tmp(lambda,j) 
        END DO
     END DO
  END DO
!********************************************************************************
  END SUBROUTINE Sector_len
!********************************************************************************
  SUBROUTINE Sector_imag( prop, eigval, eigvec, t)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:,:)                      :: prop
  REAL(idp), DIMENSION(:)                        :: eigval
  REAL(idp), DIMENSION(:,:)                      :: eigvec
  REAL(idp), DIMENSION(2,2)                      :: tmp
  REAL(idp)                                      :: t
  REAL(idp)                                      :: sum
  INTEGER                                        :: i
  INTEGER                                        :: j
  INTEGER                                        :: lambda
  prop(:,:) = 0.d0
  DO i = 1, 2
     DO lambda = 1, 2  
        tmp(i,lambda) = eigvec(i,lambda) * exp(-eigval(lambda)*t)
     END DO
  END DO
  DO i = 1, 2
     DO j = 1, 2
        DO lambda = 1 ,2
           prop(i,j) = prop(i,j) + eigvec(i,lambda) * tmp(j,lambda) 
        END DO
     END DO
  END DO
!********************************************************************************
  END SUBROUTINE Sector_imag
!********************************************************************************
!********************************************************************************
  SUBROUTINE V_Scale_z( sol, v)
  IMPLICIT NONE
  COMPLEX(idp), DIMENSION(:)                     :: sol
  REAL(idp),    DIMENSION(:)                     :: v

  sol(:) = exp(-eye * v(:) * t_b(1) ) *sol(:)
!********************************************************************************
  END SUBROUTINE V_Scale_z
!********************************************************************************
!********************************************************************************
  SUBROUTINE V_Scale_d( sol, v )
  IMPLICIT NONE
  REAL(idp),    DIMENSION(:)                     :: sol
  REAL(idp),    DIMENSION(:)                     :: v
  REAL(idp)                                      :: t
  sol(:) = exp(- v(:) * t_b(1) ) *sol(:)
!********************************************************************************
  END SUBROUTINE V_Scale_d
!********************************************************************************
  SUBROUTINE V_Scale_d_mat( sol, v )
  IMPLICIT NONE
  REAL(idp),    DIMENSION(:,:)                   :: sol
  REAL(idp),    DIMENSION(:)                     :: v
  REAL(idp)                                      :: t
  INTEGER                                        :: i
  DO i = 1, M_eig
     sol(:,i) = exp(- v(:) * t_b(1) ) *sol(:,i)
  END DO
!********************************************************************************
  END SUBROUTINE V_Scale_d_mat
!********************************************************************************
!********************************************************************************
  SUBROUTINE Exp_Off_Diag_z( so, mat, so_z, pr_z, v )
  IMPLICIT NONE
  TYPE(SPLIT_OPERATOR)                           :: so
  TYPE(MATRICES), DIMENSION(:), ALLOCATABLE      :: mat 
  TYPE(COMPLEX_H)                                :: so_z
  TYPE(COMPLEX_PROP)                             :: pr_z     
  COMPLEX(idp), DIMENSION(:)                     :: v
  INTEGER                                        :: first
  INTEGER                                        :: i
  DO trips = 1, n_trips
     first = locate(trips)
     DO i = begin(trips), locate(2), 2
        Call Mat_2_V_2(so%mat(i)%pr_z%prop,v(first:),so_z%scr(first:))
        first = first + 2 
     END DO
  END DO
!********************************************************************************
  END SUBROUTINE Exp_Off_Diag_z
!********************************************************************************
!********************************************************************************
  SUBROUTINE Exp_Off_Diag_d( so, mat, so_d, pr_d, v )
  IMPLICIT NONE
  TYPE(SPLIT_OPERATOR)                           :: so
  TYPE(MATRICES), DIMENSION(:), ALLOCATABLE      :: mat 
  TYPE(REAL_H)                                   :: so_d
  TYPE(REAL_PROP)                                :: pr_d
  TYPE(REAL_MATRIX)                              :: type_real_matrix      
  REAL(idp), DIMENSION(:)                        :: v
  INTEGER                                        :: first
  INTEGER                                        :: i
  DO trips = 1, n_trips
     first = locate(trips)
     DO i = begin(trips), mat_size, 2
        Call Mat_2_V_2(so%mat(i)%pr_d%prop,v(first:),so_d%scr(first:))
        first = first + 2 
     END DO
  END DO
!********************************************************************************
  END SUBROUTINE Exp_Off_Diag_d
!********************************************************************************
  SUBROUTINE Exp_Off_Diag_d_mat( so, mat, so_d, pr_d, v )
  IMPLICIT NONE
  TYPE(SPLIT_OPERATOR)                           :: so
  TYPE(MATRICES), DIMENSION(:), ALLOCATABLE      :: mat 
  TYPE(REAL_H)                                   :: so_d
  TYPE(REAL_PROP)                                :: pr_d
  TYPE(REAL_MATRIX)                              :: type_real_matrix      
  REAL(idp), DIMENSION(:,:)                      :: v
  INTEGER                                        :: first
  INTEGER                                        :: i
  DO trips = 1, n_trips
     first = locate(trips)
     DO i = begin(trips), mat_size, 2
        Call Mat_2_V_2(so%mat(i)%pr_d%prop,v(first:,:),so_d%scr(first:))
!        write(iout,*) 'sol after',so_d%sol
        first = first + 2 
     END DO
  END DO
!********************************************************************************
  END SUBROUTINE Exp_Off_Diag_d_mat
!********************************************************************************
!********************************************************************************
  SUBROUTINE Mat_2_V_2_z( mat, v, scr)
  IMPLICIT NONE
  Complex(idp), DIMENSION(:,:)                        :: mat
  Complex(idp), DIMENSION(:)                          :: v
  Complex(idp), DIMENSION(:)                          :: scr
  INTEGER                                             :: i
  INTEGER                                             :: j
  scr(1:2) = 0.d0
  DO i=1,2
     DO j= 1,2
        scr(j) = scr(j) + mat(j,i) * v(i)
     END DO
  END DO
  v(1:2) = scr(1:2)
!********************************************************************************
  END SUBROUTINE Mat_2_V_2_z
!********************************************************************************
  SUBROUTINE Mat_2_V_2_d( mat, v, scr)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:,:)                           :: mat
  REAL(idp), DIMENSION(:)                             :: v
  REAL(idp), DIMENSION(:)                             :: scr
  INTEGER                                             :: i
  INTEGER                                             :: j
  scr(1:2) = 0.d0
  DO i=1,2
     DO j= 1,2
        scr(j) = scr(j) + mat(j,i) * v(i)
     END DO
  END DO
  v(1:2) = scr(1:2)
!********************************************************************************
  END SUBROUTINE Mat_2_V_2_d
!********************************************************************************
!********************************************************************************
  SUBROUTINE Mat_2_Mat_2_d( mat, v, scr)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:,:)                           :: mat
  REAL(idp), DIMENSION(:,:)                           :: v
  REAL(idp), DIMENSION(:)                             :: scr
  INTEGER                                             :: i
  INTEGER                                             :: j
  INTEGER                                             :: k
  DO k = 1, M_Eig
     scr(1:2) = 0.d0
     DO i=1,2
        DO j= 1,2
           scr(j) = scr(j) + mat(j,i) * v(i,k)
        END DO
     END DO
     v(1:2,k) = scr(1:2)
  END DO
!********************************************************************************
  END SUBROUTINE Mat_2_Mat_2_d
!********************************************************************************
!********************************************************************************
  SUBROUTINE Orthonormalize_d( v_new, v_old, new, old, dim)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:,:)                           :: v_new
  REAL(idp), DIMENSION(:,:), OPTIONAL                 :: v_old
  REAL(idp)                                           :: sdot
  REAL(idp)                                           :: overlap
  INTEGER                                             :: new
  INTEGER, OPTIONAL                                   :: old
  INTEGER                                             :: dim
  INTEGER                                             :: i
  INTEGER                                             :: j
  INTEGER                                             :: k
!
! Step one is to orthogonalize the new vectors to the old ones
!
  IF ( PRESENT(v_old) ) THEN
       DO i=1,new
          DO j = 1, old
             overlap = sdot(dim,v_old(1,j),1,v_new(1,i),1)
             DO k = 1, dim
                v_new(k,i) = v_new(k,i) - v_old(k,j) * overlap
             END DO
          END DO
          overlap = 1.d0 / sqrt ( sdot(dim,v_new(1,i),1,v_new(1,i),1) )
          DO k = 1, dim
             v_new(k,i) = overlap * v_new(k,i)
          END DO
       END DO
  END IF
!!
! Now orthonormalize the new vectors to each other
!
  DO i = 1, new
     DO j = 1, i
        overlap = sdot(dim,v_new(1,j),1,v_new(1,i),1)
        DO k = 1, dim
           v_new(k,i) = v_new(k,i) - overlap * v_new(k,j)
        END DO
     END DO
     overlap = 1.d0 / sqrt ( sdot(dim,v_new(1,i),1,v_new(1,i),1) )
     DO k = 1, dim
        v_new(k,i) = overlap * v_new(k,i)
     END DO
  END DO
!********************************************************************************
  END SUBROUTINE Orthonormalize_d
!********************************************************************************
!********************************************************************************
  SUBROUTINE Orthonormalize_z( v_new, v_old, new, old, dim)
  IMPLICIT NONE
  COMPLEX(idp), DIMENSION(:,:)                        :: v_new
  COMPLEX(idp), DIMENSION(:,:), OPTIONAL              :: v_old
  COMPLEX(idp)                                        :: cdotc
  COMPLEX(idp)                                        :: overlap
  INTEGER                                             :: new
  INTEGER, OPTIONAL                                   :: old
  INTEGER                                             :: dim
  INTEGER                                             :: i
  INTEGER                                             :: j
  INTEGER                                             :: k
!
! Step one is to orthogonalize the new vectors to the old ones
!
  IF ( PRESENT(v_old) ) THEN
     DO i=1,new
        DO j = 1, old
           overlap = cdotc(dim,v_old(1,j),1,v_new(1,i),1)
           DO k = 1, dim
              v_new(k,i) = v_new(k,i) - v_old(k,j) * overlap
           END DO
        END DO
        overlap = 1.d0 / sqrt ( cdotc(dim,v_new(1,i),1,v_new(1,i),1))
        DO k = 1, dim
           v_new(k,i) = overlap * v_new(k,i)
        END DO
     END DO
  END IF
!
! Now orthonormalize the new vectors to each other
!
  DO i = 1, new
     DO j = 1, i
        overlap = cdotc(dim,v_new(1,j),1,v_new(1,i),1 )
        DO k = 1, dim
           v_new(k,i) = v_new(k,i) - overlap * v_new(k,j)
        END DO
     END DO
     overlap = 1.d0 / sqrt ( cdotc(dim,v_new(1,i),1,v_new(1,i),1))
     DO k = 1, dim
        v_new(k,i) = overlap * v_new(k,i)
     END DO
  END DO
!********************************************************************************
  END SUBROUTINE Orthonormalize_z
!********************************************************************************
!********************************************************************************
  SUBROUTINE Normalize_d( v, dim)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:)                             :: v
  REAL(idp)                                           :: sdot
  REAL(idp)                                           :: overlap
  INTEGER                                             :: dim
  overlap = 1.d0 / sqrt ( sdot(dim,v,1,v,1) )
  v = overlap * v
!********************************************************************************
  END SUBROUTINE Normalize_d
!********************************************************************************
!********************************************************************************
  SUBROUTINE Normalize_z( v, dim )
  IMPLICIT NONE
  COMPLEX(idp), DIMENSION(:)                          :: v
  COMPLEX(idp)                                        :: cdotc
  COMPLEX(idp)                                        :: overlap
  INTEGER                                             :: dim
  overlap = 1.d0 / sqrt ( cdotc(dim,v,1,v,1) )
  v = overlap * v
!********************************************************************************
  END SUBROUTINE Normalize_z
!********************************************************************************
  END MODULE SO_Module
!********************************************************************************
