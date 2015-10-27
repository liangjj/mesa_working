!***********************************************************************
! DVR_H_0_Module
!**begin prologue     DVR_H_0_Module
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Calculate the FEDVR one body matrix elements
!***description       The raw matrix KE elements are used to assemble
!***                  full one-body Kinetic energy and/or Hamiltonian
!***                  sector matrices.  Thus angular momentum terms
!***                  and any one-body potentials are added in where required.
!***                  
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      DVR_H_0_Module
!***********************************************************************
!***********************************************************************
                           MODULE DVR_H_0_Module
                           USE Data_Module
                           USE FEDVR_Shared
                           USE FEDVR_Derived_Types
                           USE Legendre
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                            INTERFACE K_0_Matrix                       
                       MODULE PROCEDURE K_0_DVR,                        &
                                        K_0_DVR_Odd,                    &
                                        K_0_DVR_Fourier,                &  
                                        K_0_DVR_Hermite,                &  
                                        K_0_DVR_Laguerre  
                            END INTERFACE K_0_Matrix
!
                            INTERFACE Add_m_Angular_Momentum                       
                       MODULE PROCEDURE Add_m_Even_Angular_Momentum,    &
                                        Add_m_Odd_Angular_Momentum
                            END INTERFACE Add_m_Angular_Momentum                       
!
                            INTERFACE Fix_BC
                       MODULE PROCEDURE Fix_BC_Nabla,                   &
                                        Fix_BC_Hamiltonian
                            END INTERFACE Fix_BC
!
REAL(idp), Private, DIMENSION(:,:), ALLOCATABLE       :: matrix
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck K_0_DVR.f
!***begin prologue     K_0_DVR
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form the basic nabla matrix elements for even FE regional matrices.
!***                   This corresponds to the l=0 angular momentum matrix.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE K_0_DVR(grid,reg_mat)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (matrices), DIMENSION(:)                    :: reg_mat
  INTEGER                                          :: ir
  CHARACTER (LEN=3)                                :: itoc
!
  DO ir = 1, nreg
     ALLOCATE(grid%reg_type_op(ir,0)%tr(1:npt(ir),1:npt(ir)))  ! Note only the zero index matrix is computed
     grid%reg_type_op(ir,0)%tr(:,:) = grid%reg_mat(ir)%ham(:,:) 
     DEALLOCATE(grid%reg_mat(ir)%ham)         
  END DO
  IF (prn(4) == .true. ) THEN
      DO ir = 1, nreg
         title = 'Basic even matrix before trim Region = '//itoc(ir)
        Call prntfmn(title,grid%reg_type_op(ir,0)%tr,                    &
                     npt(ir),npt(ir),npt(ir),npt(ir),iout,'e')
      END DO
  END IF
END SUBROUTINE K_0_DVR
!***********************************************************************
!***********************************************************************
!deck K_0_DVR_Odd.f
!***begin prologue     K_0_DVR_Odd
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form basic nabla matrix elements for odd FE regional matrices.
!***                   This corresponds to the l=1 angular momentum matrix.
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE K_0_DVR_Odd(grid,reg_mat_odd)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (odd_matrices), DIMENSION(:)                :: reg_mat_odd
  INTEGER                                          :: ir
  CHARACTER (LEN=3)                                :: itoc
!
  DO ir = 1, nreg
     ALLOCATE(grid%reg_type_op(ir,1)%tr(1:npt(ir),1:npt(ir))) ! Note only the one index matrix is computed
     grid%reg_type_op(ir,1)%tr(:,:) = grid%reg_mat_odd(ir)%ham(:,:) 
     DEALLOCATE(grid%reg_mat_odd(ir)%ham)         
  END DO
  IF (prn(4) == .true. ) THEN
      DO ir = 1, nreg
         title = 'Basic odd matrix before trim Region = '//itoc(ir)
         Call prntfmn(title,grid%reg_type_op(ir,1)%tr,                     &
                      npt(ir),npt(ir),npt(ir),npt(ir),iout,'e')
      END DO
  END IF
END SUBROUTINE K_0_DVR_Odd
!***********************************************************************
!***********************************************************************
!deck K_0_DVR_Fourier.f
!***begin prologue     K_0_DVR_Fourier
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form basic nabla matrix elements for Fourier DVR.
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE K_0_DVR_Fourier(grid,reg_mat_fourier)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (fourier_matrices), DIMENSION(:)            :: reg_mat_fourier
  INTEGER                                          :: ir
  CHARACTER (LEN=3)                                :: itoc
!
  DO ir = 1, nreg
     ALLOCATE( grid%reg_type_op(ir,0)%tr(1:npt(ir),1:npt(ir)))
     grid%reg_type_op(ir,0)%tr(:,:) = grid%reg_mat_fourier(ir)%ham(:,:) 
     DEALLOCATE(grid%reg_mat_fourier(ir)%ham) 
  END DO
  IF (prn(4) == .true. ) THEN
      DO ir = 1, nreg
         title = 'Basic Fourier matrix before trim Region = '//itoc(ir)
         Call prntfmn(title,grid%reg_type_op(ir,0)%tr,                          &
                      npt(ir),npt(ir),npt(ir),npt(ir),iout,'e')
      END DO
  END IF 
END SUBROUTINE K_0_DVR_Fourier
!***********************************************************************
!***********************************************************************
!deck K_0_DVR_Hermite.f
!***begin prologue     K_0_DVR_Hermite
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form basic nabla matrix elements for Hermiten FE regional matrices.
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE K_0_DVR_Hermite(grid,reg_mat_hermite)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (hermite_matrices), DIMENSION(:)            :: reg_mat_hermite
  INTEGER                                          :: ir
  CHARACTER (LEN=3)                                :: itoc
!
  DO ir = 1, nreg
     ALLOCATE(grid%reg_type_op(ir,0)%tr(npt(ir),npt(ir)))
     grid%reg_type_op(ir,0)%tr(:,:) = grid%reg_mat_hermite(ir)%ham(:,:) 
     DEALLOCATE(grid%reg_mat_hermite(ir)%tr) 
  END DO
  IF (prn(4) == .true. ) THEN
      DO ir = 1, nreg
         title = 'Basic Hermite matrix before trim Region = '//itoc(ir) 
         Call prntfmn(title,grid%reg_type_op(ir,0)%tr,                  &
                      npt(ir),npt(ir),npt(ir),npt(ir),iout,'e')
      END DO
  END IF
END SUBROUTINE K_0_DVR_Hermite
!***********************************************************************
!***********************************************************************
!deck K_0_DVR_Laguerre.f
!***begin prologue     K_0_DVR_Laguerre
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form basic nabla matrix elements for Laguerre regional matrices.
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE K_0_DVR_Laguerre(grid,reg_mat_laguerre)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (laguerre_matrices), DIMENSION(:)           :: reg_mat_laguerre
  INTEGER                                          :: ir
  CHARACTER (LEN=3)                                :: itoc
!
  DO ir = 1, nreg
     ALLOCATE(grid%reg_type_op(ir,0)%tr(npt(ir),npt(ir)))
     grid%reg_type_op(ir,0)%tr(:,:) =  grid%reg_mat_laguerre(ir)%ham(:,:) 
     DEALLOCATE(grid%reg_mat_laguerre(1)%ham) 
  END DO 
  IF (prn(4) == .true. ) THEN
      DO ir = 1, nreg
         title = 'Basic Laguerre matrix before trim Region = '//itoc(ir)
         Call prntfmn(title,grid%reg_type_op(ir,0)%tr,                  &
                      npt(ir),npt(ir),npt(ir),npt(ir),iout,'e')
      END DO
  END IF
END SUBROUTINE K_0_DVR_Laguerre
!***********************************************************************
!***********************************************************************
!deck Add_Potential
!***begin prologue    Add_Potential
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Input the nabla matrix, copy into the hamiltonian matrix
!***                   and add the potential.  The angular terms have already been added.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Add_Potential(grid)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  INTEGER                                          :: lm
  INTEGER                                          :: ir
  INTEGER                                          :: i
  CHARACTER (LEN=3)                                :: itoc
!
  DO lm = 0, lm_max
     DO ir = 1, nreg
        ALLOCATE(grid%reg_type_op(ir,lm)%ham(1:npt(ir),1:npt(ir)))         
        grid%reg_type_op(ir,lm)%ham(:,:) = grid%reg_type_op(ir,lm)%tr(:,:) 
        DO i = 1, npt(ir)
           grid%reg_type_op(ir,lm)%ham(i,i) = grid%reg_type_op(ir,lm)%ham(i,i)      &
                                                      +                             &
                                              reg_pot(ir)%vec(i)                    &
                                                      *                             &
                                              grid%reg_pt_wt(ir)%qr_fac(i)         
        END DO 
        grid%reg_type_op(ir,lm)%ham(:,:) = dscale * grid%reg_type_op(ir,lm)%ham(:,:) 
     END DO
  END DO
  IF (prn(4) == .true. ) THEN
      DO lm = 0, lm_max
         DO ir = 1, nreg
            title = 'Hamiltonian matrix before trim LM= '//itoc(lm)//               &
                    ' Region = '//itoc(ir)
            Call prntfmn(title,grid%reg_type_op(ir,lm)%ham,npt(ir),npt(ir),         &
                                                           npt(ir),npt(ir),iout,'e')
         END DO
      END DO
  END IF
END SUBROUTINE Add_Potential
!***********************************************************************
!***********************************************************************
!deck Add_l_Angular_Momentum
!***begin prologue     Add_l_Angular_Momentum
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Add in angular momentum to nabla matrix
!***                   The matrix with L=0 zero already is formed.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Add_l_Angular_Momentum

  SUBROUTINE Add_l_Angular_Momentum(grid)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  INTEGER                                          :: i
  INTEGER                                          :: ir
  INTEGER                                          :: l
  INTEGER                                          :: l_num
!
  DO l = 1, lm_max
     l_num = l * ( l + 1 )
     DO ir = 1, nreg
        ALLOCATE( grid%reg_type_op(ir,l)%tr(1:npt(ir),1:npt(ir)) )
        grid%reg_type_op(ir,l)%tr(:,:) = grid%reg_type_op(ir,0)%tr(:,:) 
        DO i = 1, npt(ir)
           grid%reg_type_op(ir,l)%tr(i,i) = grid%reg_type_op(ir,l)%tr(i,i)               &
                                       -                                                 &
           l_num / ( grid%reg_pt_wt(ir)%qr(i) * grid%reg_pt_wt(ir)%qr(i) )               &
                                       *                                                 &
                           grid%reg_pt_wt(ir)%qr_fac(i)
        END DO 
     END DO
  END DO
END SUBROUTINE Add_l_Angular_Momentum
!***********************************************************************
!***********************************************************************
!deck Add_m_Even_Angular_Momentum
!***begin prologue     Add_m_Even_Angular_Momentum
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Add in the m angular momentum to the nabla.
!***                   The matrix with m=zero is already formed.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Add_m_Even_Angular_Momentum
  SUBROUTINE Add_m_Even_Angular_Momentum(grid, type_even, m_start, m_skip)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (even)                                      :: type_even
  INTEGER                                          :: m
  INTEGER                                          :: m_num
  INTEGER                                          :: m_start
  INTEGER                                          :: m_skip
  INTEGER                                          :: ir
  INTEGER                                          :: i
!
  type_even%kind = 'even m nabla'
  DO m = m_start, m_max, m_skip
     m_num = m * m
     DO ir = 1, nreg
        ALLOCATE( grid%reg_type_op(ir,m)%tr(1:npt(ir),1:npt(ir)) )
        grid%reg_type_op(ir,m)%tr(:,:) = grid%reg_type_op(ir,0)%tr(:,:)   
        DO i = 1, npt(ir)
           grid%reg_type_op(ir,m)%tr(i,i) = grid%reg_type_op(ir,m)%tr(i,i)                   &
                                                      -                                      &
                                            m_num * grid%reg_pt_wt(ir)%inv_qr_fac(i)
        END DO 
     END DO
  END DO
END SUBROUTINE Add_m_Even_Angular_Momentum
!***********************************************************************
!***********************************************************************
!deck Add_m_Odd_Angular_Momentum
!***begin prologue     Add_m_Odd_Angular_Momentum
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Add in the m angular momentum to the nabla.
!***                   The matrix with m=one is already formed.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Add_m_Odd_Angular_Momentum
  SUBROUTINE Add_m_Odd_Angular_Momentum(grid, type_odd, m_start, m_skip)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (odd)                                       :: type_odd
  INTEGER                                          :: ir
  INTEGER                                          :: i
  INTEGER                                          :: m
  INTEGER                                          :: m_num
  INTEGER                                          :: m_start
  INTEGER                                          :: m_skip
!
  type_odd%kind = 'odd m nabla'
  DO m = m_start, m_max, m_skip
     m_num = m * m
     DO ir = 1, nreg
        ALLOCATE( grid%reg_type_op(ir,m)%tr(1:npt(ir),1:npt(ir)) )
        grid%reg_type_op(ir,m)%tr(:,:) = grid%reg_type_op(ir,1)%tr(:,:)   
        DO i = 1, npt(ir)
           grid%reg_type_op(ir,m)%tr(i,i) = grid%reg_type_op(ir,m)%tr(i,i)                   &
                                                      -                                      &
                                            m_num * grid%reg_pt_wt(ir)%inv_qr_fac(i)
        END DO 
     END DO
  END DO
  DO ir = 1, nreg
     DO i = 1, npt(ir)
        grid%reg_type_op(ir,1)%tr(i,i) = grid%reg_type_op(ir,1)%tr(i,i)                      &
                                          -                                                  &
                                         grid%reg_pt_wt(ir)%inv_qr_fac(i)
     END DO 
  END DO
END SUBROUTINE Add_m_Odd_Angular_Momentum
!***********************************************************************
!***********************************************************************
!deck Fix_BC_Hamiltonian
!***begin prologue     Fix_BC_Hamiltonian
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Fixing boundary conditions by removing certain basis functions
!***                   from either the first or last finite element.
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Fix_BC_Hamiltonian

  SUBROUTINE Fix_BC_Hamiltonian(grid, type_hamiltonian)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (ham)                                       :: type_hamiltonian
  INTEGER                                          :: lm
  INTEGER                                          :: n_max
  CHARACTER (LEN=3)                                :: itoc
!
  n_max = max(int_zero,npt(1),npt(nreg))
  ALLOCATE(matrix(1:n_max,1:n_max))
  IF ( drop(1) == .true.) THEN
       n_tmp(1) = npt(1)-1
       DO lm = int_zero, lm_max
          matrix(1:n_tmp(1),1:n_tmp(1)) = grid%reg_type_op(1,lm)%ham(2:npt(1),2:npt(1))
          DEALLOCATE(grid%reg_type_op(1,lm)%ham)
          ALLOCATE(grid%reg_type_op(1,lm)%ham(1:n_tmp(1),1:n_tmp(1)))
          grid%reg_type_op(1,lm)%ham(:,:) = matrix(:,:)
       END DO
       IF (prn(4) == .true. ) THEN
           DO lm = 0, lm_max
              title='First region Hamiltonian matrix LM = '//itoc(lm)//' with proper boundary conditions'
             Call prntfmn(title,grid%reg_type_op(1,lm)%ham,n_tmp(1),n_tmp(1),            &
                                                           n_tmp(1),n_tmp(1),iout,'e')
           END DO
       END IF
  END IF 
  IF ( drop(2) == .true.) THEN
       n_tmp(nreg) = npt(nreg)-1
       DO lm = int_zero, lm_max
           matrix(1:n_tmp(nreg),1:n_tmp(nreg)) = grid%reg_type_op(nreg,lm)%ham(1:npt(nreg)-1,1:npt(nreg)-1)
           DEALLOCATE(grid%reg_type_op(nreg,lm)%ham)
           ALLOCATE(grid%reg_type_op(nreg,lm)%ham(1:n_tmp(nreg),1:n_tmp(nreg)))
           grid%reg_type_op(nreg,lm)%ham(:,:) = matrix(:,:)
       END DO
       IF (prn(4) == .true. ) THEN
           DO lm = 0, lm_max
              title='Last region Hamiltonian matrix LM = '//itoc(lm)//' with proper boundary conditions'
              Call prntfmn(title,grid%reg_type_op(1,lm)%ham,n_tmp(nreg),n_tmp(nreg),         &
                                                            n_tmp(nreg),n_tmp(nreg),iout,'e')
           END DO
       END IF
  END IF 
  DEALLOCATE(matrix)
END SUBROUTINE Fix_BC_Hamiltonian
!***********************************************************************
!***********************************************************************
!deck Fix_BC_Nabla
!***begin prologue     Fix_BC_Nabla
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Fixing boundary conditions by removing certain basis functions
!***                   from either the first or last finite element.
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Fix_BC_Nabla

  SUBROUTINE Fix_BC_Nabla(grid, type_nabla)
  IMPLICIT NONE
  TYPE (coordinates)                               :: grid
  TYPE (nabla)                                     :: type_nabla
  INTEGER                                          :: n_max
  INTEGER                                          :: lm
  CHARACTER (LEN=3)                                :: itoc
!
  n_max = max(int_zero,npt(1),npt(nreg))
  ALLOCATE(matrix(1:n_max,1:n_max))
  IF ( drop(1) == .true.) THEN
       n_tmp(1) = npt(1)-1
       DO lm = int_zero, lm_max
          matrix(1:n_tmp(1),1:n_tmp(1)) = grid%reg_type_op(1,lm)%tr(2:npt(1),2:npt(1))
          DEALLOCATE(grid%reg_type_op(1,lm)%tr)
          ALLOCATE(grid%reg_type_op(1,lm)%tr(1:n_tmp(1),1:n_tmp(1)))
          grid%reg_type_op(1,lm)%tr(:,:) = matrix(:,:)
       END DO
       IF (prn(4) == .true. ) THEN
           DO lm = 0, lm_max
              title='First region nabla matrix LM = '//itoc(lm)//' with proper boundary conditions'
             Call prntfmn(title,grid%reg_type_op(1,lm)%tr,n_tmp(1),n_tmp(1),            &
                                                          n_tmp(1),n_tmp(1),iout,'e')
           END DO
       END IF
  END IF 
  IF ( drop(2) == .true.) THEN
        n_tmp(nreg) = npt(nreg)-1
       DO lm = int_zero, lm_max
           matrix(1:n_tmp(nreg),1:n_tmp(nreg)) = grid%reg_type_op(nreg,lm)%tr(1:npt(nreg)-1,1:npt(nreg)-1)
           DEALLOCATE(grid%reg_type_op(nreg,lm)%tr)
           ALLOCATE(grid%reg_type_op(nreg,lm)%tr(1:n_tmp(nreg),1:n_tmp(nreg)))
           grid%reg_type_op(nreg,lm)%tr(:,:) = matrix(:,:)
       END DO
       IF (prn(4) == .true. ) THEN
           DO lm = 0, lm_max
              title='Last region nabla matrix LM = '//itoc(lm)//' with proper boundary conditions'
              Call prntfmn(title,grid%reg_type_op(1,lm)%tr,n_tmp(nreg),n_tmp(nreg),            &
                                                           n_tmp(nreg),n_tmp(nreg),iout,'e')
           END DO
       END IF
  END IF 
  DEALLOCATE(matrix)
END SUBROUTINE Fix_BC_Nabla
!***********************************************************************
!***********************************************************************
!deck PE_DVR_Matrix.f
!***begin prologue     PE_DVR_Matrix
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate the local potential matrix for a number of
!***                   fairly standard potentials.
!***                   
!***                   
!***                   
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       PE_DVR_Matrix
  SUBROUTINE PE_DVR_Matrix(grid)
  IMPLICIT NONE
  TYPE (coordinates)                           :: grid
  INTEGER                                      :: i
  INTEGER                                      :: intkey
  CHARACTER (LEN=3)                            :: itoc
  REAL(idp)                                    :: fpkey
!
!
  DO i = 1, nreg
     ALLOCATE( reg_pot(i)%vec(npt(i)) )
     reg_pot(i)%vec(:) = zero     
!
     IF(reg_pot(i)%type == 'none') then
!
        call none(reg_pot(i)%vec,                                   &
                  grid%reg_pt_wt(i)%qr,                             &
                  npt(i),                                           &
                  prnt)
!
     ELSE IF(reg_pot(i)%type == 'well') then
!
        reg_pot(i)%well%depth=fpkey(card,'well_depth',zero,' ')
        call vwell(reg_pot(i)%vec,                                  &
                   reg_pot(i)%well%depth,                           &
                   npt(i),                                          &
                   prnt)
!
     ELSE IF(reg_pot(i)%type == 'exponential') then
!
        reg_pot(i)%exponential%amp(1)=fpkey(card,'amplitude',-1.d0,' ')
        reg_pot(i)%exponential%expnt(1)=fpkey(card,'exponent',reg_pot(i)%exponential%expnt,' ')
        call vexp(reg_pot(i)%vec,                                   &
                  grid%reg_pt_wt(i)%qr,                             &
                  reg_pot(i)%exponential%amp,                       &
                  reg_pot(i)%exponential%expnt,                     &
                  npt(i),                                           &
                  prnt)
!
     ELSE IF(reg_pot(i)%type == 'yukawa') then
!
        reg_pot(i)%yukawa%amp(1)=fpkey(card,'amplitude',-1.d0,' ')
        reg_pot(i)%yukawa%expnt(1)=fpkey(card,'exponent',reg_pot(i)%exponential%expnt,' ')
        call vyukawa(reg_pot(i)%vec,                                &
                     grid%reg_pt_wt(i)%qr,                          &
                     reg_pot(i)%yukawa%amp,                         &
                     reg_pot(i)%yukawa%expnt,                       &
                     npt(i),                                        &
                     prnt)
!
     ELSE IF(reg_pot(i)%type == 'power_exponential') then
!
        reg_pot(i)%power_exp%amp(1)=fpkey(card,'amplitude',1.d0,' ')
        reg_pot(i)%power_exp%expnt(1)=fpkey(card,'exponent',reg_pot(i)%power_exp%expnt,' ')
        reg_pot(i)%power_exp%n_p=intkey(card,'power',0,' ')
        call v_pow_exp(reg_pot(i)%vec,                               &
                       grid%reg_pt_wt(i)%qr,                         &
                       reg_pot(i)%power_exp%amp,                     &
                       reg_pot(i)%power_exp%expnt,                   &
                       reg_pot(i)%power_exp%n_p,                     &
                       npt(i),                                       &
                       prnt)
!
     ELSE IF(reg_pot(i)%type == 'sum_exponential') then
!
        call fparr(card,'amplitudes',reg_pot(i)%sum_exp%amp,2,' ')
        call fparr(card,'exponents',reg_pot(i)%sum_exp%expnt,2,' ')
        call vexp_sum(reg_pot(i)%vec,                                &
                      grid%reg_pt_wt(i)%qr,                          &
                      reg_pot(i)%sum_exp%amp,                        &
                      reg_pot(i)%sum_exp%expnt,                      &
                      npt(i),                                        &
                      prnt)
!
     ELSE IF(reg_pot(i)%type == 'coulomb') then
!
        reg_pot(i)%coulomb%charge=fpkey(card,'charge',-one,' ')
        call vcoul(reg_pot(i)%vec,                                   &
                   grid%reg_pt_wt(i)%qr,                             &
                   reg_pot(i)%coulomb%charge,                        &
                   npt(i),                                           &
                   prnt)
!
     ELSE IF(reg_pot(i)%type == 'eberlonium') then
!
        reg_pot(i)%eberlonium%charge=fpkey(card,'charge',-one,' ')
        reg_pot(i)%eberlonium%n_p=intkey(card,'power',0,' ')
        reg_pot(i)%eberlonium%amp(1)=fpkey(card,'a',1.d0,' ')
        reg_pot(i)%eberlonium%amp(2)=fpkey(card,'b',1.d0,' ')
        call v_eberlonium(reg_pot(i)%vec,                            &
                          grid%reg_pt_wt(i)%qr,                      &
                          reg_pot(i)%eberlonium%charge,              &
                          reg_pot(i)%eberlonium%amp(1),              &
                          reg_pot(i)%eberlonium%amp(2),              &
                          reg_pot(i)%eberlonium%n_p,                 &
                          npt(i),                                    &
                          prnt)
!
     ELSE IF(reg_pot(i)%type == 'inverse_r4') then
!
        call vir4(reg_pot(i)%vec,                                    &
                  grid%reg_pt_wt(i)%qr,                              &
                  npt(i),                                            &
                  prnt)
!
     ELSE IF(reg_pot(i)%type == 'rounded_well') then
!
        reg_pot(i)%well%nwell=intkey(card,'n_well',ten,' ')
        reg_pot(i)%well%awell=fpkey(card,'a_well',14.d0,' ')
        call vrwell(reg_pot(i)%vec,                                  &
                    grid%reg_pt_wt(i)%qr,                            &
                    reg_pot(i)%well%awell,                           &
                    reg_pot(i)%well%nwell,                           &
                    npt(i),                                          &
                    prnt)
!
     ELSE IF(reg_pot(i)%type == 'harmonic_oscillator') then
!
!        enter the mass and frequency in atomic units.
!         
        reg_pot(i)%harmonic_oscillator%mass=fpkey(card,'mass',1.d0,' ')
        reg_pot(i)%harmonic_oscillator%omega=fpkey(card,'omega',1.d0,' ')
        factor= reg_pot(i)%harmonic_oscillator%mass * reg_pot(i)%harmonic_oscillator%omega &
                                                    * reg_pot(i)%harmonic_oscillator%omega &
                                                    * half
        hbar=one
        write(iout,1) mass, reg_pot(i)%harmonic_oscillator%omega
        call vhmo(reg_pot(i)%vec,                                    &
                  grid%reg_pt_wt(i)%qr,                              &
                  factor,                                            &
                  npt(i),                                            &
                  prnt)
!
     ELSE IF(reg_pot(i)%type == 'anharmonic_oscillator') then
!
        call vanhmo(reg_pot(i)%vec,                                  &
                    grid%reg_pt_wt(i)%qr,                            &
                    npt(i),                                          &
                    prnt)
!
     ELSE IF(reg_pot(i)%type == 'expres') then
!
        call fparr(card,'amplitude',reg_pot(i)%expres%amp,2,' ')
        call fparr(card,'exponent',reg_pot(i)%expres%expnt,2,' ')
        reg_pot(i)%expres%shift=fpkey(card,'exponent_shift',zero,' ')
        call vres(reg_pot(i)%vec,                                    &
                  grid%reg_pt_wt(i)%qr,                              &
                  reg_pot(i)%expres%amp,                             &
                  reg_pot(i)%expres%expnt,                           &
                  reg_pot(i)%expres%shift,                           &
                  npt(i),                                            &
                  prnt)
!
     ELSE IF(reg_pot(i)%type == 'periodic') then
!
        reg_pot(i)%periodic%n_scale=intkey(card,'n_i',10,' ')         
        reg_pot(i)%periodic%e_c=fpkey(card,'e_c',.001d0,' ')         
        reg_pot(i)%periodic%awell=reg_pot(i)%periodic%n_scale/reg_pot(i)%periodic%e_c
        call vperiod(reg_pot(i)%vec,                                 &
                     grid%reg_pt_wt(i)%qr,                           &
                     reg_pot(i)%periodic%awell,                      &
                     npt(i),                                         &
                     prnt)
!
     ELSE IF(reg_pot(i)%type == 'spheroidal') then
!
!    In the spheroidal case, this is not actually the potential but a scaled version which
!    takes into account the volume factor in the integral which is non-separable in the two
!    coordinates.
!
        IF ( grid%label == 'eta') THEN
             reg_pot(i)%vec(:) = R_ab * ( Z_b - Z_a ) * grid%reg_pt_wt(i)%qr(:) * grid%reg_pt_wt(i)%inv_qr_fac(:)
        ELSE IF ( grid%label == 'xi' ) THEN
             reg_pot(i)%vec(:) = R_AB * ( Z_a + Z_b ) * grid%reg_pt_wt(i)%qr(:) * grid%reg_pt_wt(i)%inv_qr_fac(:)
        END IF
     ELSE
        Call lnkerr('screwed up potential. quit')
     END IF
!
  END DO
  IF(prn(3)) then
     DO i = 1, nreg
        write(iout,2) reg_pot(i)%type
        title='potential matrix elements for region = '//itoc(i)
        call prntrm(title,reg_pot(i)%vec,npt(i),1,npt(i),1,iout)
     END DO
  END IF
!
!
  dscale = - half * hbar * hbar / mass
  IF(units == 'atomic_units') then
     dscale = - half
  END IF     
!
  IF(reg_pot(i)%type == 'spheroidal') then
     dscale = - R_ab * quarter
  END IF
!
 1 Format(/,1x,'oscillator mass      = ',e15.8, &
          /,1x,'oscillator-frequency = ',e15.8)
 2 FORMAT(/,20x,'Potential Type = ',a32)
!
END SUBROUTINE PE_DVR_Matrix
!***********************************************************************
!***********************************************************************
           END MODULE DVR_H_0_Module
!***********************************************************************
!***********************************************************************
