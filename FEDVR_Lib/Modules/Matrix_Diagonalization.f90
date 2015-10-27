!***********************************************************************
! Matrix_Diagonalization
!**begin prologue     Matrix_Diagonalization
!**date written       090219   (yymmdd)
!**revision date               (yymmdd)
!**keywords           DVR, FEDVR
!**
!**author             schneider, b. i.(nsf)
!**source             DVR Library
!**purpose            Diagonalize the final DVR matrix
!***                  
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      Matrix_Diagonalization
!***********************************************************************
!***********************************************************************
                           MODULE Matrix_Diagonalization
                           USE DVR_H_0_Module
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Diagonalize_Global_Matrices
!***begin prologue     Diagonalize_Global_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Diagonalize the Hamiltonian
!***                   
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Diagonalize_Global_Matrices
!
  SUBROUTINE Diagonalize_Global_Matrices(grid)
  IMPLICIT NONE
  TYPE (coordinates)                             :: grid
  INTEGER                                        :: count
  INTEGER                                        :: i
  INTEGER                                        :: j
  INTEGER                                        :: lm
  INTEGER                                        :: info
  CHARACTER(LEN=3)                               :: itoc
  LOGICAL                                        :: dollar
  LOGICAL                                        :: logkey
  INTEGER                                        :: intkey
  CHARACTER(LEN=80)                              :: chrkey
!  
  write(iout,*)
  write(iout,*) '                              Diagonalizing the Full Matrix'
  write(iout,*)
  ALLOCATE(dvr_mat(0)%eigenvectors(1:physical_points,1:physical_points),                &
           dvr_mat(0)%eigenvalues(1:physical_points),                                   &
           dvr_mat(0)%work(3*physical_points),                                          &
           dvr_mat(0)%lower( physical_points*(physical_points+1)/2 ) )
  DO lm = 0, size
!
!    Use the packed form of the diagonalization routine
!
     count = 0
     DO i = 1, physical_points
        DO j = 1, i
           count = count + 1
           dvr_mat(0)%lower(count) = dvr_mat(lm)%ham(i,j)
        END DO
     END DO
     DEALLOCATE(dvr_mat(lm)%ham)
!
     Call dspev('v','u',physical_points,dvr_mat(0)%lower,              &
                dvr_mat(0)%eigenvalues,dvr_mat(0)%eigenvectors,        &
                physical_points,dvr_mat(0)%work,info)
!
!
     title='eigenvalues of Hamiltonian for angular quantum number = '//itoc(lm)
     CALL prntrm(title,dvr_mat(0)%eigenvalues,physical_points,1,physical_points,1,iout)
     IF(prn(10)) THEN
        title='eigenvectors'    
        CALL prntrm(title,dvr_mat(0)%eigenvectors,physical_points,physical_points,                       &
                                                  physical_points,physical_points,iout)
     END IF
     write(iout,*) 'Writing eigenvalues and eigenvectors to disk'
     Call IOsys('write real "'//title(1:len)//' eigenvalues" to '//FEDVR_File,            &
                 physical_points,dvr_mat(0)%eigenvalues,0,' ')
     Call IOsys('write real "'//title(1:len)//' eigenvectors" to '//FEDVR_File,           &
                 physical_points*physical_points,dvr_mat(0)%eigenvectors,0,' ')
  END DO
  DEALLOCATE(dvr_mat(0)%eigenvectors,                                                   &
             dvr_mat(0)%eigenvalues,                                                    &
             dvr_mat(0)%work,                                                           &
             dvr_mat(0)%lower )
END SUBROUTINE Diagonalize_Global_Matrices
!***********************************************************************
!***********************************************************************
           END MODULE Matrix_Diagonalization
!***********************************************************************
!***********************************************************************
