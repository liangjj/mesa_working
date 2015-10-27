  MODULE Crank_Nicholson_Module
  USE Data
  USE Derived_Types
  USE Matrix_Print
  USE Time_Propagation_Subroutines_Module
  IMPLICIT NONE
!
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                            INTERFACE CN_Propagation                       
                       MODULE PROCEDURE CN_Length_Propagation,                           &
                                        CN_Velocity_Propagation
                            END INTERFACE CN_Propagation
!____________________________________________________________________________________________!
!____________________________________________________________________________________________!
!
                                  Contains
!********************************************************************************
!********************************************************************************
  SUBROUTINE CN_Driver 
  IMPLICIT NONE
  TYPE(CRANK_NICHOLSON)                             :: cn
  TYPE(CN_LENGTH)                                   :: cn_len
  TYPE(CN_VELOCITY)                                 :: cn_vel
  INTEGER                                           :: len
!
  IF (Method == 'cn_length') THEN
        Call CN_Propagation(cn,cn_len)                       
  ELSE IF (Method == 'cn_velocity') THEN
      Call CN_Propagation(cn,cn_vel)     
  ELSE IF (Method == 'cn_length_and_cn_velocity') THEN    
      Call CN_Propagation(cn,cn_len)    
      Call CN_Propagation(cn,cn_vel)                                     
  ELSE
      Call lnkerr('bad call to crank-nicholson')
  END IF
!********************************************************************************
  END SUBROUTINE CN_Driver 
!********************************************************************************
  SUBROUTINE CN_Length_Propagation (cn,cn_len)
  IMPLICIT NONE
  TYPE(Complex_Vector)                       :: type_complex_vector
  TYPE(CRANK_NICHOLSON)                      :: cn
  TYPE(CN_LENGTH)                            :: cn_len
  TYPE(ENE_HAM)                              :: hamiltonian
  INTEGER                                    :: i
  INTEGER                                    :: j
!
  ALLOCATE (cn%ham_l(1:m_size), cn%ham_u(1:m_size), cn%ham_d(1:m_size),                    &
            cn%sol(1:m_size), cn%lhs_d(1:m_size), cn%lhs_u(1:m_size), cn%lhs_l(1:m_size),  &
            cn%rhs_d(1:m_size), cn%rhs_u(1:m_size), cn%rhs_l(1:m_size),                    &
            work_d(1:m_size), V_t(1:m_size) )
!
  write(iout,*)
  cn_len%length_title='Starting the CN Length Propagation Routine'
  Write(iout,*) cn_len%length_title
  pulse_time  =  delta_t * m_size * 0.5d0
  write(iout,*) '         Pulse Time = ',pulse_time
! 
  cn%sol(1:m_size)    = Z (1:m_size,1)    
  cn%ham_l(1:m_size)  = E (1:m_size)
  cn%ham_u(1: m_size) = E (1:m_size)
!  
  call cpu_time (t_i)
  time = 0.0d0
  E_flag = .false.
  DO i = 1, NO_Time_Steps
     t_l = time + .5d0*delta_t
     IF ( E_flag /= .true. ) THEN
          Call Comp_E (t_l)
          Call TD_Pot(t_l)
          cn%ham_d(1:m_size) = D(1:m_size) + V_t(1:m_size)
     ELSE
          cn%ham_d(1:m_size) = D(1:m_size) + V_t(1:m_size)
     END IF
     time  = time + delta_t
     CALL CN_Solver(cn)
     IF (Prnt(2) == .true. ) Then
         write(iout,1) i
         Write(iout,2) ( x(j), real(cn%sol(j)), imag(cn%sol(j)), j=1,m_size)
     END IF
     Call Get_Probabilities(cn%sol,work_d)
  END DO
  write(iout,*) '         Final Time = ',time
  call cpu_time (t_f)
  write (iout,*) '***** Length Solution *****'
  IF ( t_f - t_i < 60 ) THEN
     write (iout,*) 'Time elapsed: ', (t_f - t_i),'seconds'
  ELSE
     write (iout,*) 'Time elapsed: ', (t_f - t_i)/real(60) , 'minutes'
  END IF
  IF (info == 0) THEN
  
!------------------------------ TEST 1 -----------------------------------!
!------------------ COMPUTING < LENGTH | H | LENGTH > --------------------!
!----IF AT E = 0, ENERGY IS EQUAL TO THE GROUND STATE ENERGY, THEN THE----!
!---------------------------TEST IS SUCCESSFULL---------------------------! 
     call Compute_Energies (hamiltonian,cn%sol,cn%ham_d)
     write (iout,*)
     write (iout,*) '"TEST 1"' 
     write (iout,*) 'Energy AT E_Field = ', E_0, ': ', Energy
!-------------------------------------------------------------------------!
!-------------------------------------------------------------------------!

     Call Print_Matrix(type_complex_vector,cn%sol,title='CN Length Gauge Solution Vector:', &
                       frmt='fr')
  ELSE 
     print *, 'info is not zero in CN Length Gauge Solution'
  END IF   
1 format(/,10x,'Time Step = ',i5,10x,'x',10x,'real part of wf', 10x 'imaginary part of wf') 
2 Format(10x,25x,d15.8,10x,d15.8,10x,d15.8)
  DEALLOCATE (cn%ham_l, cn%ham_u, cn%ham_d, cn%sol, cn%lhs_d, cn%lhs_u, cn%lhs_l,           &
              cn%rhs_d, cn%rhs_u, cn%rhs_l )
!********************************************************************************
  END SUBROUTINE CN_Length_Propagation 
!********************************************************************************
  SUBROUTINE CN_Velocity_Propagation (cn, cn_vel)
  IMPLICIT NONE
  TYPE(Complex_Vector)                          :: type_complex_vector
  TYPE(CRANK_NICHOLSON)                         :: cn
  TYPE(CN_VELOCITY)                             :: cn_vel
  TYPE(ENE_HAM)                                 :: hamiltonian
  INTEGER                                       :: i
  INTEGER                                       :: j
!
  write(iout,*)
  cn_vel%velocity_title='Starting the CN Velocity Propagation Routine'
  Write(iout,*) cn_vel%velocity_title
  pulse_time  =  delta_t * m_size * 0.5d0
  write(iout,*) '          Pulse Time = ',pulse_time
  quad_size = 4
  ALLOCATE (cn%ham_l(1:m_size), cn%ham_u(1:m_size), cn%ham_d(1:m_size),                    &
            cn%sol(1:m_size), cn%lhs_d(1:m_size), cn%lhs_u(1:m_size), cn%lhs_l(1:m_size),  &
            cn%rhs_d(1:m_size), cn%rhs_u(1:m_size), cn%rhs_l(1:m_size), pt(1:quad_size),   &
            wt(1:quad_size), work_d(1:quad_size), V_t(1:m_size) )

  call cpu_time (t_i)
  time = 0.0d0
  t_l = time 
  t_u = t_l + .5d0*delta_t
  A=0.d0
  CALL Comp_A(A,t_l,t_u,pt,wt,work_d)
  cn%sol(:)  = Z(:,1)
  E_flag = .false.
  DO i =1, No_time_steps
     IF (E_flag /= .true. ) THEN
         cn%ham_d (1:m_size) = D(1:m_size) + A*A*0.5d0
         cn%ham_l(1:m_size)  = E(1:m_size) - eye * A * 0.5d0 / Step_Size
         cn%ham_u(1:m_size)  = E(1:m_size) + eye * A * 0.5d0 / Step_Size
         A_old = A
         t_l = t_u
         t_u = t_l + delta_t
         CALL Comp_A (A, t_l, t_u,  pt, wt, work_d)
     ELSE
         cn%ham_d (1:m_size) = D(1:m_size) 
         cn%ham_l(1:m_size)  = E(1:m_size) 
         cn%ham_u(1:m_size)  = E(1:m_size) 
         time = time + delta_t
     ENDIF
     CALL CN_Solver(cn)
     IF (Prnt(2) == .true. ) Then
         write(iout,1) i
         Write(iout,2) ( x(j), real(cn%sol(j)), imag(cn%sol(j)), j=1,m_size)
     END IF
     time = time + delta_t
  END DO
  write(iout,*) '          Final Time = ',time
  call cpu_time (t_f)
  write (iout,*) '***** Velocity Solution *****'
  IF ( t_f - t_i < 60 ) THEN
     write (iout,*) 'Time elapsed: ', (t_f - t_i),'seconds'
  ELSE
     write (iout,*) 'Time elapsed: ', (t_f - t_i)/real(60),'minutes'
  END IF

  IF (info == 0) THEN
!------------------------------ TEST 1 -----------------------------------!
!---------------- COMPUTING < VELOCITY | H | VELOCITY > ------------------!
!----IF AT E = 0, ENERGY IS EQUAL TO THE GROUND STATE ENERGY, THEN THE----!
!---------------------------TEST IS SUCCESSFULL---------------------------!
!      cn%sol=exp(-eye*A_old*x) * cn%sol
!      call Test_Energy (cn%sol,cn%ham_d)
      cn%sol=exp(-eye*A_old*x) * cn%sol
      call Compute_Energies (hamiltonian,cn%sol,cn%ham_d)
      write (iout,*)
      write (iout,*) '"TEST 1"'
      write (iout,*) 'Energy AT E_field = ', E_0, ': ', Energy
!-------------------------------------------------------------------------!

     Call Print_Matrix(type_complex_vector,cn%sol,title='CN Velocity Gauge Solution Vector:', &
                       frmt='fr')
  ELSE 
      write(iout,*) 'info is not zero in Velocity Gauge Solution'
  END IF 
1 format(/,10x,'Time Step = ',i5,10x,'x',10x,'real part of wf', 10x 'imaginary part of wf') 
2 Format(10x,25x,d15.8,10x,d15.8,10x,d15.8)
  DEALLOCATE (cn%ham_l, cn%ham_u, cn%ham_d, cn%sol, cn%lhs_d, cn%lhs_u, cn%lhs_l,             &
              cn%rhs_d, cn%rhs_u, cn%rhs_l, pt, wt, work_d )  
!********************************************************************************
  END SUBROUTINE  CN_Velocity_Propagation
!********************************************************************************
!******************************************************************************** 
!                            CRANK-NICOLSON
!********************************************************************************
  SUBROUTINE CN_Solver(cn)
  IMPLICIT NONE
  TYPE(CRANK_NICHOLSON)                      :: cn
  INTEGER                                    :: i
  COMPLEX(idp)                               :: t_1
  COMPLEX(idp)                               :: t_2
!
  t_1  = eye * delta_t * 0.5d0
  DO i = 1, m_size
     t_2           = cn%ham_d(i) * t_1
     cn%lhs_d(i) = 1.d0 + t_2
     cn%rhs_d(i) = 1.d0 - t_2
     cn%lhs_u(i) = t_1 * cn%ham_u(i)
     cn%lhs_l(i) = t_1 * cn%ham_l(i)
     cn%rhs_u(i) = -t_1 * cn%ham_u(i)
     cn%rhs_l(i) = -t_1 * cn%ham_l(i)
  END DO
!
  cn%ham_d(1)      =      cn%rhs_d(1) * cn%sol(1)               &
                                      +                         &
                          cn%rhs_u(1) * cn%sol(2)
  cn%ham_d(m_size) =      cn%rhs_l(m_size-1) * cn%sol(m_size-1) &
                                           +                    &
                        cn%rhs_d(m_size)   * cn%sol(m_size)
   DO i = 2, m_size-1
      cn%ham_d(i) = cn%rhs_l(i-1) * cn%sol(i-1)                 &
                                +                               &
                  cn%rhs_d(i)   * cn%sol(i)                     &
                                +                               &
                  cn%rhs_u(i+1) * cn%sol(i+1)
  END DO
  cn%sol(1:m_size) = cn%ham_d(1:m_size)
!..........................................................................!
  CALL ZGTSV( m_size, 1, cn%lhs_l, cn%lhs_d, cn%lhs_u, cn%sol, m_size, info)
!*******************************************************************************            
  END SUBROUTINE CN_Solver
!*******************************************************************************
!********************************************************************************
  END MODULE Crank_Nicholson_Module
!********************************************************************************
