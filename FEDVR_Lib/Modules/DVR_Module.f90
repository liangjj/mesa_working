!***********************************************************************
! DVR_Module
!**begin prologue     DVR_Module
!**date written       090219   (yymmdd)
!**revision date               (yymmdd)
!**keywords           DVR, FEDVR
!**
!**author             schneider, b. i.(nsf)
!**source             DVR Library
!**purpose            This is the driver module to calculate the FEDVR functions 
!***                  and matrix elements.  The actual work is done in the 
!***                  subroutines that are called.
!***                  in a FEDVR basis
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      DVR_Module
!***********************************************************************
!***********************************************************************
                           MODULE DVR_Module
                           USE Data_Module
                           USE FEDVR_Shared
                           USE FEDVR_Derived_Types
                           USE DVR_Polynomials_Module
                           USE DVR_Kinetic_Energy_Module
                           USE DVR_H_0_Module
                     INTEGER                      :: len_1
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck KE_DVR_Matrices
!***begin prologue     KE_DVR_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Driving routine to compute the raw sector matrix elements
!***                   of the kinetic energy operator.  These are done for each sector
!***                   with no regard for boundary conditions or normalization.
!***                   For certain coordinates even and odd m must be treated differently.
!***                   The real work is done in the DVR_Kinetic_Energy_Module where the
!***                   quadrature is performed by Kinetic_Energy.
!***  
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE KE_DVR_Matrices(grid)
  IMPLICIT NONE
  TYPE (coordinates)             :: grid
!
  len=lenth(grid%label)
  IF (keyword == 'cartesian') THEN
!
      IF(typwt == 'one' .or. typwt == 'legendre') THEN      
         ALLOCATE(grid%reg_mat(1:nreg))
         Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
      ELSE IF (typwt == 'hermite') THEN
         ALLOCATE(grid%reg_mat_hermite(1:nreg))
         Call Kinetic_Energy(grid, grid%reg_mat_hermite, grid%reg_poly)
      END IF
!
  ELSE IF(keyword =='spherical') THEN
!
      IF(grid%label(1:len) == 'r') THEN
         IF(typwt == 'one' .or. typwt == 'legendre') THEN
            ALLOCATE(grid%reg_mat(1:nreg))
            Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
         ELSE IF (typwt == 'laguerre') THEN
            ALLOCATE(grid%reg_mat_laguerre(1))
            Call Kinetic_Energy(grid, grid%reg_mat_laguerre, grid%reg_poly)
         ELSE IF (typwt == 'spherical_hermite') THEN        
            ALLOCATE(grid%reg_mat_hermite(1:nreg))
            Call Kinetic_Energy(grid, grid%reg_mat_hermite, grid%reg_poly)
         END IF
!
      ELSE IF( grid%label(1:len) == 'theta') THEN
!
!        Compute the Even KE.
!
         pre_factor = - one
         ALLOCATE(grid%reg_mat(1:nreg))
         Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
!
!            Compute the Odd KE.
!
         IF (m_max > 0 ) THEN
             ALLOCATE(grid%reg_mat_odd(1:nreg))
             Call Kinetic_Energy(grid, grid%reg_mat_odd, grid%reg_poly)
         END IF
      END IF
!
  ELSE IF(keyword =='cylindrical') THEN
!
      IF(grid%label(1:len) == 'rho') THEN
         ALLOCATE(grid%reg_mat(1:nreg))
         Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
      END IF
!
  ELSE IF (keyword == 'spheroidal') THEN
!
      ALLOCATE(grid%reg_mat(1:nreg))
      IF (grid%label(1:len) == 'eta') THEN
          pre_factor = - one
      ELSE IF(grid%label(1:len) == 'xi') THEN
          pre_factor = one
      END IF
!
!     Compute the even kinetic energy matrix elements.
!
      Call Kinetic_Energy(grid, grid%reg_mat, grid%reg_poly)
      IF ( m_max > 0 ) THEN
           ALLOCATE(grid%reg_mat_odd(1:nreg))
!
!      Compute the odd kinetic energy matrix elements.
!
           Call Kinetic_Energy(grid, grid%reg_mat_odd, grid%reg_poly)
      END IF
  ELSE IF (keyword == 'fourier') THEN
!
      ALLOCATE( grid%reg_mat_fourier(1))
      Call Kinetic_Energy(grid, grid%reg_mat_fourier)
!
!
  END IF
!
!
END SUBROUTINE KE_DVR_Matrices
!***********************************************************************
!***********************************************************************
!deck Final_Functions
!***begin prologue     Final_Functions
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Take the basic interpolation polynomials in each sector
!***                   and renormalize them.  The renormalization factors are computed
!***                   so that the polyomials at beginning and end of the sector, the
!***                   bridge functions, are properly normalized. See the normaliztion
!***                   routine for details. 
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE Final_Functions(grid)
  IMPLICIT NONE
  TYPE (coordinates)             :: grid
  INTEGER                        :: ir
  INTEGER                        :: i
!
!
!                  Normalize the functions.
!
  DO ir = 1, nreg 
     write(iout,1) ir
     DO i = 1, npt(ir)     
        grid%reg_poly(ir)%pr(:,i)                                                    &
                                  =                                                  &
        grid%reg_poly(ir)%pr(:,i) * grid%reg_poly(ir)%normalization(i)
        grid%reg_poly(ir)%dpr(:,i)                                                   &
                                  =                                                  &
        grid%reg_poly(ir)%dpr(:,i) * grid%reg_poly(ir)%normalization(i)
        grid%reg_poly(ir)%ddpr(:,i)                                                  &
                                  =                                                  &
        grid%reg_poly(ir)%ddpr(:,i) * grid%reg_poly(ir)%normalization(i)
     END DO
     IF (prn(3) == .true. ) THEN
         title = 'normalized polynomials'
         Call prntfmn(title,grid%reg_poly(ir)%pr,npt(ir),npt(ir),                    &
                                                 npt(ir),npt(ir),iout,'e')
         title = 'first derivative of normalized polynomials'
         Call prntfmn(title,grid%reg_poly(ir)%dpr,npt(ir),npt(ir),                   &
                                                  npt(ir),npt(ir),iout,'e')
         title = 'second derivative of normalized polynomials'
         Call prntfmn(title,grid%reg_poly(ir)%ddpr,npt(ir),npt(ir),                  &
                                                   npt(ir),npt(ir),iout,'e')
     END IF
  END DO
!
1 Format(/,10x,'Region = ',i4)
END SUBROUTINE Final_Functions
!***********************************************************************
!***********************************************************************
!deck Final_KE_DVR_Matrices
!***begin prologue     Final_KE_DVR_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute the renormalized sector matrix elements
!***                   of the kinetic energy operator.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Final_KE_DVR_Matrices
!
  SUBROUTINE Final_KE_DVR_Matrices(grid)
  IMPLICIT NONE
  TYPE (coordinates)             :: grid
!
!
  len=lenth(grid%label)
!
  IF (keyword == 'cartesian') THEN
!
      IF(typwt == 'one' .or. typwt == 'legendre') THEN      
         Call Renormalization(grid,grid%reg_mat)
      ELSE IF (typwt == 'hermite') THEN
         Call Renormalization(grid,grid%reg_mat_hermite)
      END IF
!
  ELSE IF(keyword =='spherical') THEN
!
      IF(grid%label(1:len) == 'r') THEN
!
         IF(typwt == 'one' .or. typwt == 'legendre') THEN
            Call Renormalization(grid,grid%reg_mat)
         ELSE IF (typwt == 'laguerre') THEN
            Call Renormalization(grid,grid%reg_mat_laguerre)
         ELSE IF (typwt == 'hermite') THEN
            Call Renormalization(grid,grid%reg_mat_hermite)
         ELSE IF (typwt == 'spherical_hermite') THEN
            Call Renormalization(grid,grid%reg_mat_hermite)
         END IF
!
      ELSE IF( grid%label(1:len) == 'theta') THEN
!
!        Compute the Even KE.
!
            Call Renormalization(grid,grid%reg_mat)
!
!        Compute the Odd polynomials.
!
         IF (m_max > 0 ) THEN
             Call Renormalization(grid,grid%reg_mat_odd)
         END IF
      END IF
!
  ELSE IF(keyword =='cylindrical') THEN
!
      IF(grid%label(1:len) == 'rho') THEN
         Call Renormalization(grid,grid%reg_mat)
      END IF
!
  ELSE IF (keyword == 'spheroidal') THEN
!
!
!     Compute the even kinetic energy matrix elements.
!
         Call Renormalization(grid,grid%reg_mat)
      IF ( m_max > 0 ) THEN
!
!      Compute the odd kinetic energy matrix elements.
!
           Call Renormalization(grid, grid%reg_mat_odd)
      END IF
  ELSE IF (keyword == 'fourier') THEN
           Call Renormalization(grid, grid%reg_mat_fourier)
  END IF
!
END SUBROUTINE Final_KE_DVR_Matrices
!***********************************************************************
!***********************************************************************
!deck H_0.f
!***begin prologue     H_0
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            The renormalized sector matrix elements of the kinetic energy 
!***                   operator with proper boundary conditions are combined with
!***                   any one body potential.  Then the angular momentum parts are
!***                   added and the matrices "trimmed" if there are boundary conditions
!***                   at the end points which need to be satisfied.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE H_0(grid)
  IMPLICIT NONE
  TYPE (coordinates)                   :: grid
  TYPE (even)                          :: type_even
  TYPE (odd)                           :: type_odd
  TYPE (nabla)                         :: type_nabla
  TYPE (ham)                           :: type_hamiltonian
  INTEGER                              :: i
  INTEGER                              :: lm
  CHARACTER (LEN=3)                    :: itoc
!
!          The first step is to take the sector matrix elements which have been computed
!          without regard to angular momentun or potential energy terms and to store them
!          as specific KE matrices for each angular momentum.  The raw elements are then
!          deallocated.  Basically that is all K_0_Matrix does
!          The second step is to add in any angular momentum terms where required.
!          The third step is to copy those matrices into the hamiltonian array, adding in the
!          potential.
!          The fourth step is take the result of step 3. and trim the first and last dvr element if
!          needed to get the physical matrices.  At the end of this routine we have what we need
!          to proceed.
!
!          Store the original size of the first and last element in these two variables.  They will be
!          reset in Fix_BC if that routine is called and then reset later.
!
  ALLOCATE(n_tmp(1:nreg))
  n_tmp(1:nreg) = npt(1:nreg)
  IF (keyword == 'cartesian') THEN
!
      ALLOCATE(grid%reg_type_op(1:nreg,0:lm_max))
      IF(typwt == 'one' .or. typwt == 'legendre') THEN      
         Call K_0_Matrix(grid, grid%reg_mat)
      ELSE IF(typwt == 'hermite') THEN
         Call K_0_Matrix(grid, grid%reg_mat_hermite)
      END IF
!   
!                No angular terms, just make the hamiltonian
!
      lm_max = int_zero
!
      IF (form_hamiltonian) THEN
          write(iout,*) 'Forming Cartesian Hamiltonian'
!
          Call Add_Potential(grid)
!
!         Trim if necesssary and replace the first and last element sizes by their trimmed values
!
          IF (drop(1) == .true. .or. drop(2) == .true. ) THEN
              Call Fix_BC(grid,type_hamiltonian)
          END IF
          IF (prn(4) == .true. ) THEN
              DO lm = 0, lm_max
                 DO i = 1, nreg
                    title = 'Cartesian Hamiltonian matrix Region '//itoc(i)//' with proper boundary conditions'
                    Call prntfmn(title,grid%reg_type_op(i,lm)%ham,n_tmp(i),n_tmp(i),           &
                                                                  n_tmp(i),n_tmp(i),iout,'e')
                 END DO
              END DO
          END IF
      END IF
!
      IF (form_nabla) THEN
          write(iout,*) 'Forming Cartesian Nabla'
!
!         Trim if necesssary and replace the first and last element sizes by their trimmed values
!
          IF (drop(1) == .true. .or. drop(2) == .true. ) THEN
              Call Fix_BC(grid,type_nabla)
          END IF
          IF (prn(4) == .true. ) THEN
              DO lm = 0, lm_max
                 DO i = 1, nreg
                    title = 'Cartesian Nabla matrix Region '//itoc(i)//' with proper boundary conditions'
                    Call prntfmn(title,grid%reg_type_op(i,lm)%tr,n_tmp(i),n_tmp(i),            &
                                                                 n_tmp(i),n_tmp(i),iout,'e')
                 END DO
              END DO
          END IF
      ELSE
          DO lm = 0, lm_max
             DO i = 1,nreg
                DEALLOCATE(grid%reg_type_op(i,lm)%tr)
             END DO
          END DO
      END IF
!
  ELSE IF(keyword =='spherical') THEN
!
      IF(grid%label(1:len) == 'r') THEN
         ALLOCATE(grid%reg_type_op(1:nreg,0:l_max))
         IF(typwt == 'one' .or. typwt == 'legendre') THEN
            Call K_0_Matrix(grid, grid%reg_mat)
         ELSE IF (typwt == 'laguerre') THEN
            Call K_0_Matrix(grid, grid%reg_mat_laguerre)
         ELSE IF (typwt == 'hermite') THEN
            Call K_0_Matrix(grid, grid%reg_mat_hermite)
         ELSE IF (typwt == 'spherical_hermite') THEN
            Call K_0_Matrix(grid, grid%reg_mat_hermite)
         END IF
         lm_max = l_max
         Call Add_l_Angular_Momentum(grid)
         IF (form_hamiltonian) THEN
             Call Add_Potential(grid)
             write(iout,*)
             write(iout,*) 'Radial Hamiltonian constructed'
             IF (drop(1) == .true. .or. drop(2) == .true. ) THEN
                 Call Fix_BC(grid,type_hamiltonian)
             END IF
             IF (prn(4) == .true. ) THEN
                 DO lm = 0, lm_max
                    DO i = 1, nreg
                       title = 'Spherical Radial Hamiltonian matrix LM  '//itoc(lm)//          &
                               ' Region '//itoc(i)//' with proper boundary conditions'
                       Call prntfmn(title,grid%reg_type_op(i,lm)%ham,n_tmp(i),n_tmp(i)    ,    &
                                                                     n_tmp(i),n_tmp(i),iout,'e')
                    END DO
                 END DO
             END IF
         END IF
!
         IF (form_nabla) THEN
             write(iout,*)
             write(iout,*) 'Radial Nabla constructed'
             IF (drop(1) == .true. .or. drop(2) == .true. ) THEN
                 Call Fix_BC(grid,type_nabla)
             END IF
             IF (prn(4) == .true. ) THEN
                 DO lm = 0, lm_max
                    DO i = 1, nreg
                       title = 'Spherical Radial Nabla matrix LM  '//itoc(lm)//                &
                               ' Region '//itoc(i)//' with proper boundary conditions'
                       Call prntfmn(title,grid%reg_type_op(i,lm)%tr,n_tmp(i),n_tmp(i),         &
                                                                    n_tmp(i),n_tmp(i),iout,'e')
                    END DO
                 END DO
             END IF
         ELSE
             DO lm = 0, lm_max
                DO i = 1, nreg
                   DEALLOCATE(grid%reg_type_op(i,lm)%tr)
                END DO
             END DO
         END IF
      ELSE IF( grid%label(1:len) == 'theta') THEN
!
!
         ALLOCATE(grid%reg_type_op(1:nreg,0:m_max))
!
!        Compute the Even KE.
!
         lm_max = m_max
         Call K_0_Matrix(grid, grid%reg_mat)
         Call Add_m_Angular_Momentum(grid, type_even, int_two, int_two)
         write(iout,*) 'Finished theta '//type_even%kind
!
!        Compute the Odd KE.
!
         IF (m_max > 0 ) THEN
!
             Call K_0_Matrix(grid, grid%reg_mat_odd)
             Call Add_m_Angular_Momentum(grid, type_odd, int_three, int_two)
             write(iout,*) 'FInished theta '//type_odd%kind
!
         END IF
         IF (form_hamiltonian) THEN
             write(iout,*) 'Forming Spherical Theta Hamiltonian'
             Call Add_Potential(grid)
             IF (drop(1) == .true. .or. drop(2) == .true. ) THEN
                 Call Fix_BC(grid,type_hamiltonian)
             END IF
             IF (prn(4) == .true. ) THEN
                 DO lm = 0, lm_max
                    DO i = 1, nreg
                       title = 'Theta Hamiltonian matrix LM = '//itoc(lm)//                    &
                               ' Region = '//itoc(i)//' with proper boundary conditions'
                       Call prntfmn(title,grid%reg_type_op(i,lm)%ham,n_tmp(i),n_tmp(i),        &
                                                                     n_tmp(i),n_tmp(i),iout,'e')
                    END DO
                 END DO
             END IF
         END IF
!
         IF (form_nabla) THEN
             write(iout,*) 'Forming Spherical Theta Nabla'
             IF (drop(1) == .true. .or. drop(2) == .true. ) THEN
                 Call Fix_BC(grid,type_nabla)
             END IF
             IF (prn(4) == .true. ) THEN
                 DO lm = 0, lm_max
                    DO i = 1, nreg
                       title = 'Theta Nabla matrix LM = '//itoc(lm)//                          &
                               ' Region = '//itoc(i)//' with proper boundary conditions'
                       Call prntfmn(title,grid%reg_type_op(i,lm)%tr,n_tmp(i),n_tmp(i),         &
                                                                    n_tmp(i),n_tmp(i),iout,'e')
                    END DO
                 END DO
             END IF
         ELSE
             DO lm = 0, lm_max
                DO i = 1, nreg
                   DEALLOCATE(grid%reg_type_op(i,lm)%tr)
                END DO
             END DO
         END IF
      END IF
!
  ELSE IF(keyword =='cylindrical') THEN
!
      IF(grid%label(1:len) == 'rho') THEN
         ALLOCATE(grid%reg_type_op(1:nreg,0:m_max))         
         Call K_0_Matrix(grid, grid%reg_mat)
         lm_max = m_max
         Call Add_m_Angular_Momentum(grid, type_even, int_one, int_one)
         write(iout,*) 'Finished cylindrical '//type_even%kind
         drop(1) = .false.
         IF (form_hamiltonian) THEN
             write(iout,*) 'Forming Cylindrical Rho Hamiltonian'
             Call Add_Potential(grid)
             IF (drop(2) == .true. ) THEN
                 Call Fix_BC(grid,type_hamiltonian)
             END IF
             IF (prn(4) == .true. ) THEN
                 DO lm = 0, lm_max
                    DO i = 1, nreg
                       title = 'Rho Hamiltonian matrix LM = '//itoc(lm)//                      &
                               ' Region = '//itoc(i)//' with proper boundary conditions'
                       Call prntfmn(title,grid%reg_type_op(i,lm)%ham,n_tmp(i),n_tmp(i),        &
                                                                     n_tmp(i),n_tmp(i),iout,'e')
                    END DO
                 END DO
             END IF
         END IF
!
         IF (form_nabla) THEN
             write(iout,*) 'Forming Cylindrical Rho Nabla'
             IF (drop(2) == .true. ) THEN
                 Call Fix_BC(grid,type_nabla)
             END IF
             IF (prn(4) == .true. ) THEN
                 DO lm = 0, lm_max
                    DO i = 1, nreg
                       title = 'Rho Nabla matrix LM = '//itoc(lm)//                            &
                               ' Region = '//itoc(i)//' with proper boundary conditions'
                       Call prntfmn(title,grid%reg_type_op(i,lm)%tr,n_tmp(i),n_tmp(i),         &
                                                                    n_tmp(i),n_tmp(i),iout,'e')
                    END DO
                 END DO
             END IF
         ELSE   
             DO lm = 0, lm_max
                DO i = 1, nreg
                   DEALLOCATE(grid%reg_type_op(i,lm)%tr)
                END DO
             END DO
         END IF
      END IF
  ELSE IF (keyword == 'spheroidal') THEN
!
      IF (grid%label(1:len) == 'eta') THEN
          pre_factor = - one
      ELSE IF(grid%label(1:len) == 'xi') THEN
          pre_factor = one
      END IF
!
!     Compute the even kinetic energy matrix elements.
!
      ALLOCATE(grid%reg_type_op(1:nreg,0:m_max))         
      Call K_0_Matrix(grid, grid%reg_mat)
      lm_max = m_max
      Call Add_m_Angular_Momentum(grid, type_even, int_two, int_two)
      write(iout,*) 'Finished spheroidal '//type_even%kind

!
!      Compute the odd kinetic energy matrix elements.
!
      IF ( m_max > 0 ) THEN
           Call K_0_Matrix(grid, grid%reg_mat_odd)
           Call Add_m_Angular_Momentum(grid, type_odd, int_three, int_two)
           write(iout,*) 'Finished spheroidal '//type_odd%kind
      END IF
      IF (form_hamiltonian) THEN
          write(iout,*) 'Forming Spheroidal '//grid%label(1:len)//' Hamiltonian'
          Call Add_Potential(grid)
          IF (drop(1) == .true. .or. drop(2) == .true. ) THEN
              Call Fix_BC(grid,type_hamiltonian)
          END IF
          IF (prn(4) == .true. ) THEN
              DO lm = 0, lm_max
                 DO i = 1, nreg
                    title = grid%label(1:len)//' Hamiltonian matrix LM = '//itoc(lm)//         &
                            ' Region = '//itoc(i)//' with proper boundary conditions'
                    Call prntfmn(title,grid%reg_type_op(i,lm)%ham,n_tmp(i),n_tmp(i),           &
                                                                  n_tmp(i),n_tmp(i),iout,'e')
                 END DO
              END DO
          END IF
      END IF
      IF (form_nabla) THEN
          write(iout,*) 'Forming Spheroidal '//grid%label(1:len)//' Nabla'
          IF (drop(1) == .true. .or. drop(2) == .true. ) THEN
              Call Fix_BC(grid,type_nabla)
          END IF
          IF (prn(4) == .true. ) THEN
              DO lm = 0, lm_max
                 DO i = 1, nreg
                    title = grid%label(1:len)//' Nabla matrix LM = '//itoc(lm)//               &
                            ' Region = '//itoc(i)//' with proper boundary conditions'
                    Call prntfmn(title,grid%reg_type_op(i,lm)%tr,n_tmp(i),n_tmp(i),            &
                                                                 n_tmp(i),n_tmp(i),iout,'e')
                 END DO
              END DO
          END IF
      ELSE   
          DO lm = 0, lm_max
             DO i = 1, nreg
                DEALLOCATE(grid%reg_type_op(i,lm)%tr)
             END DO
          END DO
      END IF
  ELSE IF (keyword == 'fourier') THEN
!
      lm_max = int_zero
      ALLOCATE(grid%reg_type_op(1:nreg,0:0))
      Call K_0_Matrix(grid, grid%reg_mat_fourier)
      IF (form_hamiltonian) THEN
          write(iout,*) 'Forming Fourier Hamiltonian'
          Call Add_Potential(grid)
          IF (prn(4) == .true. ) THEN
              DO lm = 0, lm_max
                 DO i = 1, nreg
                    title = 'Fourier  Hamiltonian matrix Region = '//itoc(i)//' with proper boundary conditions'
                    Call prntfmn(title,grid%reg_type_op(i,0)%ham,n_tmp(i),n_tmp(i),            &
                                                                 n_tmp(i),n_tmp(i),iout,'e')
                 END DO
              END DO
          END IF
      END IF
      IF (form_nabla) THEN
          write(iout,*) 'Forming Fourier Nabla'
          IF (prn(4) == .true. ) THEN
              DO i = 1, nreg
                 title = 'Fourier  Hamiltonian matrix Region = '//itoc(i)//' with proper boundary conditions'
                 Call prntfmn(title,grid%reg_type_op(i,0)%tr,n_tmp(i),n_tmp(i),                &
                                                             n_tmp(i),n_tmp(i),iout,'e')
              END DO
          END IF
      ELSE      
          DO i = 1, nreg
             DEALLOCATE(grid%reg_type_op(i,0)%tr)
          END DO
      END IF
  END IF
  npt(1:nreg) = n_tmp(1:nreg)
  DEALLOCATE(n_tmp)
!
  DO i = 1, nreg
     DEALLOCATE(reg_pot(i)%vec)
  END DO
  DEALLOCATE(reg_pot)
!
END SUBROUTINE H_0
!***********************************************************************
!***********************************************************************
           END MODULE DVR_Module
!***********************************************************************
!***********************************************************************
