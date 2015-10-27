!***********************************************************************
! Matrix_Scale_and_Assemble
!**begin prologue     Matrix_Scale_and_Assemble
!**date written       090219   (yymmdd)
!**revision date               (yymmdd)
!**keywords           DVR, FEDVR
!**
!**author             schneider, b. i.(nsf)
!**source             DVR Library
!**purpose            Assemble the final DVR matrix
!***                  in a FEDVR basis
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      Matrix_Scale_and_Assemble
!***********************************************************************
!***********************************************************************
                           MODULE Matrix_Scale_and_Assemble
                           USE DVR_H_0_Module
                           USE Data_Module
  INTEGER :: len_1
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Form_Matrix
!***begin prologue     Form_Matrix
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Assemble the full Hamiltonian matrix elements from the sector matrices
!***description        The routine assembles the full matrix from the sector matrices after transforming 
!***                   them to standard form.  The final result is the global matrix ready for further use.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Form_Matrix
!
  SUBROUTINE Form_Matrix(grid)
  IMPLICIT NONE
  TYPE (coordinates)                             :: grid
  INTEGER                                        :: i
  CHARACTER(LEN=3)                               :: itoc
  LOGICAL                                        :: dollar
  LOGICAL                                        :: logkey
  INTEGER                                        :: intkey
  CHARACTER(LEN=80)                              :: chrkey
!
!         Some coordinate systems and coordinates have a metric which needs to scale the Hamiltonian
!
  write(iout,*)
  write(iout,*) '                              Forming the Final Full Matrix'
  write(iout,*)
  len=lenth(grid%label)
  size=max(0,l_max,m_max) 
  ALLOCATE(dvr_mat(0:size))
  IF ( keyword == 'spherical' ) THEN
       IF ( grid%label(1:len) == 'r') THEN
            Call Transform_Matrix_to_Standard_Form(grid)          
       END IF
  END IF
  IF ( keyword == 'cylindrical' ) THEN
       IF ( grid%label(1:len) == 'rho') THEN
            Call Transform_Matrix_to_Standard_Form(grid)          
       END IF
  END IF
!
  Call Form_Global_Matrix(grid,grid%label(1:len))
!
!
END SUBROUTINE Form_Matrix
!***********************************************************************
!***********************************************************************
!deck Transform_Matrix_to_Standard_Form
!***begin prologue     Transform_Matrix_to_Standard_Form
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            In certain coordinate systems there is a simple metric
!                      due to the structure of the kinetic energy matrix
!                      which turn a standard eigenvalue problem to what looks
!                      like a general eigenvalue problem.  This is trivially
!                      removed by a simple diagonal pre and post multiplication.
!                      This unitary transformation also makes the definition of
!                      the vectors different and this needs to be accoubnted for later.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE Transform_Matrix_to_Standard_Form(grid)
  IMPLICIT NONE
  TYPE (coordinates)             :: grid
  INTEGER                        :: lm
  INTEGER                        :: i
  INTEGER                        :: j
  CHARACTER (LEN=3)              :: itoc
!
  IF (form_hamiltonian) THEN
      write(iout,*)
      write(iout,*) 'Forming Scaled Hamiltonian'
      DO lm = 0, lm_max
         DO i = 1, nreg
            DO j = 1, npt(i)
!
!              Pre and post multiplication
!
               grid%reg_type_op(i,lm)%ham(j,1:j) =                                       &
                                                   grid%reg_pt_wt(i)%inv_sqrt_qr_fac(j)  &
                                                           *                             &
                                                   grid%reg_type_op(i,lm)%ham(j,1:j)     &
                                                           *                             &
                                                   grid%reg_pt_wt(i)%inv_sqrt_qr_fac(1:j)
               grid%reg_type_op(i,lm)%ham(1:j,j) = grid%reg_type_op(i,lm)%ham(j,1:j) 
            END DO
         END DO
      END DO
      IF (prn(4) == .true. ) THEN
          DO lm = 0, lm_max
             DO i = 1, nreg
                title = 'Scaled Hamiltonian matrix LM = '//itoc(lm)//' Region '//itoc(i)
                Call prntfmn(title,grid%reg_type_op(i,lm)%ham,npt(i),npt(i),           &
                                                              npt(i),npt(i),iout,'e')
             END DO
          END DO
      END IF
  END IF
!
  IF (form_nabla) THEN
      write(iout,*)
      write(iout,*) 'Forming Scaled Nabla'
      DO lm = 0, lm_max
         DO i = 1, nreg
            DO j = 1, npt(i)
!
!              Pre and post multiplication
!
               grid%reg_type_op(i,lm)%tr(j,1:j) =                                        &
                                                   grid%reg_pt_wt(i)%inv_sqrt_qr_fac(j)  &
                                                           *                             &
                                                   grid%reg_type_op(i,lm)%tr(j,1:j)      &
                                                           *                             &
                                                   grid%reg_pt_wt(i)%inv_sqrt_qr_fac(1:j)
               grid%reg_type_op(i,lm)%tr(1:j,j) = grid%reg_type_op(i,lm)%tr(j,1:j) 
            END DO
         END DO
      END DO
      IF (prn(4) == .true. ) THEN
          DO lm = 0, lm_max
             DO i = 1, nreg
                title = 'Scaled Nabla matrix LM = '//itoc(lm)//' Region '//itoc(i)
                Call prntfmn(title,grid%reg_type_op(i,lm)%tr,npt(i),npt(i),              &
                                                              npt(i),npt(i),iout,'e')
             END DO
          END DO
      END IF
  END IF
!
END SUBROUTINE Transform_Matrix_to_Standard_Form
!***********************************************************************
!***********************************************************************
!deck Form_Global_Matrix.f
!***begin prologue     Form_Global_Matrix
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form the global matrices.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE Form_Global_Matrix(grid,label)
  IMPLICIT NONE
  TYPE (coordinates)                      :: grid
  CHARACTER(LEN=*)                        :: label
  INTEGER                                 :: i
  INTEGER                                 :: lm
  INTEGER                                 :: first
  INTEGER                                 :: last
  CHARACTER(LEN=3)                        :: itoc
!  
  IF (form_hamiltonian) THEN
      DO lm = 0, lm_max
         ALLOCATE(dvr_mat(lm)%ham(1:physical_points,1:physical_points))
         dvr_mat(lm)%ham(:,:) = zero
         first = 1
         DO i = 1, nreg
            last = first + npt(i) - 1
            dvr_mat(lm)%ham(first:last,first:last) = grid%reg_type_op(i,lm)%ham(1:npt(i),1:npt(i))         
            first = last
            DEALLOCATE(grid%reg_type_op(i,lm)%ham)
         END DO    
         IF (prn(9) == .true. ) THEN
             title = 'Global Hamiltonian matrix L or M = '//itoc(lm)
             Call prntfmn(title,dvr_mat(lm)%ham,physical_points,physical_points,           &
                                                physical_points,physical_points,iout,'e')
         END IF
         len=lenth(keyword)
         len_1=lenth(label)
         title=keyword(1:len)//'_'//label(1:len_1)//'_'//itoc(lm)
         len=lenth(title)
         write(iout,*) 'Writing Global Hamiltonian = ',title(1:len)//' to disk'
         matrix_title='"1_electron_'//itoc(lm)//' hamiltonian matrix file title"'
         len_1=lenth(matrix_title)
         Call IOsys('write character '//matrix_title(1:len_1)//' to '//FEDVR_File,0,0,0,title(1:len))
         Call IOsys('write real "'//title(1:len)//' hamiltonian" to '//FEDVR_File,           &
                     physical_points*physical_points,dvr_mat(lm)%ham,0,' ')
      END DO
  END IF
  IF (form_nabla) THEN
      DO lm = 0, lm_max
         ALLOCATE(dvr_mat(lm)%tr(1:physical_points,1:physical_points))
         dvr_mat(lm)%tr(:,:) = zero
         first = 1
         DO i = 1, nreg
            last = first + npt(i) - 1
            dvr_mat(lm)%tr(first:last,first:last) = grid%reg_type_op(i,lm)%tr(1:npt(i),1:npt(i))         
            first = last
            DEALLOCATE(grid%reg_type_op(i,lm)%tr)
         END DO    
         IF (prn(9) == .true. ) THEN
             title = 'Global Nabla matrix L or M = '//itoc(lm)
             Call prntfmn(title,dvr_mat(lm)%tr,physical_points,physical_points,           &
                                               physical_points,physical_points,iout,'e')
         END IF
         len=lenth(keyword)
         len_1=lenth(label)
         title=keyword(1:len)//'_'//label(1:len_1)//'_'//itoc(lm)
         len=lenth(title)
         write(iout,*) 'Writing Global Nabla = ',title(1:len)//' to disk'
         matrix_title='"1_electron_'//itoc(lm)//' nabla matrix file title"'
         len_1=lenth(matrix_title)
         Call IOsys('write character '//matrix_title(1:len_1)//' to '//FEDVR_File,0,0,0,title(1:len))
         Call IOsys('write real "'//title(1:len)//' nabla" to '//FEDVR_File,           &
                     physical_points*physical_points,dvr_mat(lm)%tr,0,' ')
      END DO
  END IF
END SUBROUTINE Form_Global_Matrix
!***********************************************************************
!***********************************************************************
           END MODULE Matrix_Scale_and_Assemble
!***********************************************************************
!***********************************************************************
