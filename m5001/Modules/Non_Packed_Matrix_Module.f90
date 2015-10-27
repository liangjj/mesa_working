!***********************************************************************
                           MODULE Non_Packed_Matrix_Module
!
                           USE Atomic_Matrices
                           USE Pack_Global
                           USE Global_Time_Propagation_Module
                           USE Pack_Hamiltonian_Module
                           USE Iterative_Global
                           USE Matrix_Utility_Module
                           USE input_output
                           USE Matrix_Print
                           USE FEDVR_Shared, ONLY : file_loc
!
                           IMPLICIT NONE
            CHARACTER(LEN=1600)                  :: data_card
            CHARACTER(LEN=80)                    :: pass_data
            CHARACTER(LEN=80)                    :: local_title
            INTEGER, DIMENSION(10)               :: len
            TYPE(REAL_MATRIX)                    :: type_real_matrix
            TYPE(REAL_VECTOR)                    :: type_real_vector
            TYPE(COMPLEX_VECTOR)                 :: type_complex_vector
            TYPE(COMPLEX_MATRIX)                 :: type_complex_matrix
!***********************************************************************
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
                          INTERFACE Small_Matrix
                    MODULE PROCEDURE Small_Matrix_d,                                      &
                                    Small_Matrix_z
                          END INTERFACE Small_Matrix
!
                          INTERFACE Input_Full_Matrices
                    MODULE PROCEDURE Input_Full_Matrices_d,                              &
                                    Input_Full_Matrices_z
                          END INTERFACE Input_Full_Matrices
!
                          INTERFACE Reformat_Matrices
                   MODULE PROCEDURE Reformat_Matrices_d,                                 &
                                    Reformat_Matrices_z
                          END INTERFACE Reformat_Matrices
!
                          INTERFACE Input_Channel_Matrices
                   MODULE PROCEDURE Input_Channel_Matrices_d,                            &
                                    Input_Channel_Matrices_z
                          END INTERFACE Input_Channel_Matrices
!
                          INTERFACE Channel_Sub_Matrix
                   MODULE PROCEDURE Channel_Sub_Matrix_d,                                &
                                    Channel_Sub_Matrix_z
                          END INTERFACE Channel_Sub_Matrix
!
                          INTERFACE Compute_Eigenstates
                   MODULE PROCEDURE Compute_Eigenstates_d,                               &
                                    Compute_Eigenstates_z                           
                          END INTERFACE Compute_Eigenstates
!
                          INTERFACE Compute_Inverse
                   MODULE PROCEDURE Compute_Inverse_d,                                   &
                                    Compute_Inverse_z                           
                          END INTERFACE Compute_Inverse
!
                          INTERFACE New_Matrix
                   MODULE PROCEDURE New_Matrix_d,                                        &
                                    New_Matrix_z
                          END INTERFACE New_Matrix
!
                          INTERFACE scale_matrix
                   MODULE PROCEDURE scale_matrix_d,                                      &
                                    scale_matrix_z
                          END INTERFACE scale_matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           Contains
!***********************************************************************
!***********************************************************************
!deck Small_Matrix_d
!***begin prologue     Small_Matrix_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Read in and store in main memory the triangle of a matrix.
!***                   Optionally one may reformat the matrices
!***                   using channel labelling and/or diagonalize the matrix.
!***
!***references
!***routines called
!***end prologue       Small_Matrix_d
!
  SUBROUTINE Small_Matrix_d(prop_mat_d)
  IMPLICIT NONE
  TYPE(REAL_PROP_VAR)                   :: prop_mat_d
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  INTEGER                               :: i
!
!!
  ALLOCATE(rowlab(n3d),collab(n3d))
  DO i=1,n3d
     rowlab(i) = 'Row '//itoc(i)
     collab(i) = 'Col '//itoc(i)
  END DO
!
!     
   Write(iout,1)
   ALLOCATE(prop_mat_d%triangle_hamiltonian(tri_size))
   IF(non_orth) THEN      
      ALLOCATE(prop_mat_d%triangle_overlap(tri_size))
   END IF
!
!     Input the matrices
!
   Call Input_Full_Matrices(prop_mat_d%triangle_hamiltonian, prop_mat_d%triangle_overlap) 
!
!     Reformat them by channel if needed.
!
   IF(reformatting_control == 'channel_format') THEN
      IF(channel_labels_only) THEN
         Write(iout,2)
         Call Reformat_Matrices(prop_mat_d%triangle_hamiltonian, prop_mat_d%triangle_overlap)
         n3d = new_size
         tri_size = new_tri_size
      ELSE 
         Write(iout,3)
         Call Input_Channel_Matrices(prop_mat_d%triangle_hamiltonian, prop_mat_d%triangle_overlap)
         n3d = new_size
         tri_size = new_tri_size
      END IF
   END IF
!
!     This section of code is just used to make explicit diagonalizations
!     for testing purposes, if needed.
!
   IF( diagonalize_only) THEN
       Write(iout,4)
       Call Compute_Eigenstates(prop_mat_d%triangle_hamiltonian, prop_mat_d%triangle_overlap) 
   ELSE IF(test_s_inverse) THEN
       Write(iout,5)
       Call Compute_Inverse(prop_mat_d%triangle_overlap) 
   END IF
!
!
  DEALLOCATE(rowlab,collab)
1  Format(/,30x,'***** 1. Input Matrices *****')
2  Format(/,30x,'***** 2. Reformat Matrices *****')
3  Format(/,30x,'***** 3. Compute Channel Matrices *****')
4  Format(/,30x,'***** 4. Eigenstate Test *****')
5  Format(/,30x,'***** 4. Inverse Test *****')
END SUBROUTINE Small_Matrix_d
!***********************************************************************
!***********************************************************************
!deck Small_Matrix_z
!***begin prologue     Small_Matrix_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Read in and store in main memory the triangle of a matrix.
!***                   Optionally one may reformat the matrices
!***                   using channel labelling and/or diagonalize the matrix.
!***
!***references
!***routines called
!***end prologue       Small_Matrix_z
!
  SUBROUTINE Small_Matrix_z(prop_mat_z)
  IMPLICIT NONE
  TYPE(COMPLEX_PROP_VAR)                :: prop_mat_z
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  INTEGER                               :: i
!
!
  ALLOCATE(rowlab(n3d),collab(n3d))
  DO i=1,n3d
     rowlab(i) = 'Row '//itoc(i)
     collab(i) = 'Col '//itoc(i)
  END DO
!
!     
   Write(iout,1)
   ALLOCATE(prop_mat_z%triangle_hamiltonian(tri_size))
   IF(non_orth) THEN      
      ALLOCATE(prop_mat_z%triangle_overlap(tri_size))
   END IF
!
!     Input the matrices
!
   Call Input_Full_Matrices(prop_mat_z%triangle_hamiltonian, prop_mat_z%triangle_overlap) 
!
!     Reformat them by channel if needed.
!
   IF(reformatting_control == 'channel_format') THEN
      IF(channel_labels_only) THEN
         Write(iout,2)
         Call Reformat_Matrices(prop_mat_z%triangle_hamiltonian, prop_mat_z%triangle_overlap)
         n3d = new_size
         tri_size = new_tri_size
      ELSE 
         Write(iout,3)
         Call Input_Channel_Matrices(prop_mat_z%triangle_hamiltonian, prop_mat_z%triangle_overlap)
         n3d = new_size
         tri_size = new_tri_size
      END IF
   END IF
!
!     This section of code is just used to make explicit diagonalizations
!     for testing purposes, if needed.
!
   IF( diagonalize_only) THEN
       Write(iout,4)
       Call Compute_Eigenstates(prop_mat_z%triangle_hamiltonian, prop_mat_z%triangle_overlap) 
   ELSE IF(test_s_inverse) THEN
       Write(iout,5)
       Call Compute_Inverse(prop_mat_z%triangle_overlap) 
   END IF
!
!
  DEALLOCATE(rowlab,collab)
1  Format(/,30x,'***** 1. Input Matrices *****')
2  Format(/,30x,'***** 2. Reformat Matrices *****')
3  Format(/,30x,'***** 3. Compute Channel Matrices *****')
4  Format(/,30x,'***** 4. Eigenstate Test *****')
5  Format(/,30x,'***** 4. Inverse Test *****')
END SUBROUTINE Small_Matrix_z
!***********************************************************************
!***********************************************************************
!Deck Input_Full_Matrices_d
!***begin prologue     Input_Full_Matrices_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Perform the actual read of the lower triangle of a real 
!***                   symmetric matrix form disk.
!***
!***references
!***routines called
!***end prologue       Input_Full_Matrices_d
!
  SUBROUTINE Input_Full_Matrices_d (ham,s)
  IMPLICIT NONE
  TYPE(REAL_MATRICES)                   :: mat_d
  REAL(idp), DIMENSION(:)               :: ham
  REAL(idp), DIMENSION(:)               :: s
  CHARACTER (LEN=8)                     :: itoc
  REAL(idp)                             :: fpkey
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: IOSTAT
!
  IF( mat%device == 'from_input') THEN
      IF(non_orth) THEN
         count = 0
         DO i=1,n3d
            call fparr(data_card,'s_'//itoc(i),s(count+1),i,' ')
            count = count + i
         END DO
         IF(print_cc) THEN          
            count = 0
            local_title='Input Overlap Matrix'
            write(iout,1) local_title
             DO i=1,n3d
                write(iout,2) i
                write(iout,3) (s(j), j=count+1, count + i )
                count = count + i
             END DO       
          END IF
          Call scale_matrix(s,smallest,largest,count,'overlap')
      END IF
!
!            Do Cholesky Decomposition and put it on Disk.
!
      ALLOCATE(mat_d%upper(tri_size))
      mat_d%upper(1:tri_size) = s(1:tri_size)
      call dpptrf('u',n3d,mat_d%upper,info)
      IF(info /= 0) THEN
         write(iout,*) 'Singular Overlap Matrix:Quit'
         Call lnkerr('Quit.  Singular Overlap Matrix')
      END IF
      write(iout,*) 'Opening File to Write Full Cholesky Decomposition'
      Call pakstr(cholesky_file_name,len(1))
      file_loc=File_Directory(4)(1:len_dir(4))//'/'//cholesky_file_name(1:len(1))
      Call pakstr(file_loc,len(2))
      Call IOsys('open cholesky_decomposition as new',0,0,0,file_loc(1:len(2)))
      Call IOsys('write real cholesky_factor to cholesky_decomposition',       &
                  tri_size,mat_d%upper,0,' ')
      write(iout,*) 'Closing File Holding Full Cholesky Decomposition'
      Call IOsys('close cholesky_decomposition',0,0,0,' ')
      count = 0
      DO i=1,n3d
         call fparr(data_card,'h_'//itoc(i),ham(count+1),i,' ')      
         count = count + i
      END DO
      Call scale_matrix(ham,smallest,largest,count,'hamiltonian')
      IF(print_cc) THEN          
         count = 0
         local_title='Input Hamiltonian Matrix'
         write(iout,1) local_title
         DO i=1,n3d
            write(iout,2) i
            write(iout,3) (ham(j), j=count+1, count + i )
            count = count + i
         END DO       
      END IF
  ELSE
      Call pakstr(mat%device,len(1))
      file_loc=File_Directory(4)(1:len_dir(4))//'/'//mat%device(1:len(1))
      Call pakstr(file_loc,len(2))
      OPEN(UNIT=50,FILE=file_loc(1:len(2)),ACCESS='sequential', FORM='unformatted',IOSTAT=IOSTAT,STATUS='old')
      READ(50)
      IF(non_orth) THEN      
         count = 0
         DO i=1,n3d
            READ(50) ( s(j), j=count+1, count + i )
            count = count + i
         END DO
         Call scale_matrix(s,smallest,largest,count,'overlap')
         IF(print_cc) THEN          
            local_title='Input Overlap Matrix'
            write(iout,1) local_title
            count = 0
            DO i=1,n3d
               write(iout,2) i
               write(iout,3) ( s(j), j=count+1, count + i )
               count = count + i
            END DO
         END IF
!
!            Do Cholesky Decomposition and put it on Disk.
!
         ALLOCATE(mat_d%upper(tri_size))
         mat_d%upper(1:tri_size) = s(1:tri_size)
         call dpptrf('u',n3d,mat_d%upper,info)
         IF(info /= 0) THEN
            write(iout,*) 'Singular Overlap Matrix:Quit'
            Call lnkerr('Quit.  Singular Overlap Matrix')
         END IF
         write(iout,*) 'Opening File to Write Full Cholesky Decomposition'
         Call pakstr(cholesky_file_name,len(1))
         file_loc=File_Directory(4)(1:len_dir(4))//'/'//'Cholesky_Decomposition'
         Call pakstr(file_loc,len(2))
         Call IOsys('open cholesky_decomposition as new',0,0,0,file_loc(1:len(2)))
         Call IOsys('write real cholesky_factor to cholesky_decomposition',       &
                     tri_size,mat_d%upper,0,' ')
         write(iout,*) 'Closing File Holding Full Cholesky Decomposition'
      END IF
      count = 0
      DO i=1,n3d
         READ(50) ( ham(j), j=count+1, count + i )      
         count = count + i
      END DO
      Call scale_matrix(ham,smallest,largest,count,'hamiltonian')
      IF(print_cc) THEN          
         local_title='Input Hamiltonian Matrix'
         write(iout,1) local_title
         count = 0
         DO i=1,n3d
            write(iout,2) i
            write(iout,3) (ham(j), j=count+1, count + i )
            count = count + i
         END DO       
      END IF
  END IF
1 Format(a80)
2 Format(/,5x,'Row = ',i4)
3 Format( (15x,5f10.5) )
!***********************************************************************
  END SUBROUTINE Input_Full_Matrices_d
!***********************************************************************
!***********************************************************************
!deck Input_Full_Matrices_z
!***begin prologue     Input_Full_Matrices_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Perform the actual read of the lower triangle of a Hermitian matrix 
!***                   from disk.
!***
!***references
!***routines called
!***end prologue       Input_Full_Matrices_z
!
  SUBROUTINE Input_Full_Matrices_z(ham,s)
  IMPLICIT NONE
  TYPE(COMPLEX_MATRICES)                :: mat_z
  COMPLEX(idp), DIMENSION(:)            :: ham 
  COMPLEX(idp), DIMENSION(:)            :: s 
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  CHARACTER (LEN=8)                     :: kewrd
  REAL(idp)                             :: fpkey
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: lenth
  INTEGER                               :: IOSTAT
  REAL(idp), DIMENSION(:), ALLOCATABLE  :: rinput
!
  ALLOCATE(rinput(n3d))
  IF( mat%device == 'from_input') THEN
      IF(non_orth) THEN     
         count = 0 
         DO i=1,n3d
            call fparr(data_card,'s_'//itoc(i),rinput(1:i),i,' ')
            s(count+1:count+i) = rinput(1:i) 
            count = count + i
         END DO
         Call scale_matrix(s,smallest,largest,count,'overlap')
         IF(print_cc) THEN          
            local_title='Input Overlap Matrix'
            write(iout,1) local_title
            count = 0
            DO i=1,n3d
               write(iout,2) i
               write(iout,3) (s(j), j=count+1, count + i )
               count = count + i
            END DO       
         END IF
         ALLOCATE(mat_z%upper(tri_size))
         mat_z%upper(1:tri_size) = s(1:tri_size)
         call zpptrf('u',n3d,mat_z%upper,info)
         IF(info /= 0) THEN
            write(iout,*) 'Singular Overlap Matrix:Quit'
            Call lnkerr('Quit.  Singular Overlap Matrix')
         END IF
         write(iout,*) 'Opening File to Write Full Cholesky Decomposition'
         Call pakstr(cholesky_file_name,len(1))
         file_loc=File_Directory(4)(1:len_dir(4))//'/'//cholesky_file_name(1:len(1))
         Call pakstr(file_loc,len(2))
         Call IOsys('open cholesky_decomposition as new',0,0,0,file_loc(1:len(2)))
         Call IOsys('write real cholesky_factor to cholesky_decomposition',       &
                     2*tri_size,mat_z%upper,0,' ')
         write(iout,*) 'Closing File Holding Full Cholesky Decomposition'
      END IF
      count = 0
      DO i=1,n3d
         call fparr(data_card,'h_'//itoc(i),rinput(1:i),i,' ')
         ham(count+1:count+i) = rinput(1:i) 
         count = count + i
      END DO
      Call scale_matrix(ham,smallest,largest,count,'hamiltonian')
      IF(print_cc) THEN          
         count = 0
         local_title='Input Hamiltonian Matrix'
         write(iout,1) local_title
         DO i=1,n3d
            write(iout,2) i
            write(iout,3) (ham(j), j=count+1, count + i )
            count = count + i
         END DO       
      END IF
  ELSE 
      Call pakstr(mat%device,len(1))
      file_loc=File_Directory(4)(1:len_dir(4))//'/'//mat%device(1:len(1))
      Call pakstr(file_loc,len(2))
      OPEN(UNIT=50,FILE=file_loc(1:len(2)),ACCESS='sequential', FORM='unformatted',IOSTAT=IOSTAT,STATUS='old')
      READ(50)
      IF(non_orth) THEN      
        count = 0
         DO i=1,n3d
            READ(50) rinput(1:i)
            s(count+1:count+i) = rinput(1:i) 
            count = count + i
         END DO
         Call scale_matrix(s,smallest,largest,count,'overlap')
         IF(print_cc) THEN          
            local_title='Input Overlap Matrix'
            write(iout,1) local_title
            count = 0
            DO i=1,n3d
               write(iout,2) i
               write(iout,3) ( s(j), j=count+1, count + i )
               count = count + i
            END DO
         END IF
         ALLOCATE(mat_z%upper(tri_size))
         mat_z%upper(1:tri_size) = s(1:tri_size)
         call zpptrf('u',n3d,mat_z%upper,info)
         IF(info /= 0) THEN
            write(iout,*) 'Singular Overlap Matrix:Quit'
            Call lnkerr('Quit.  Singular Overlap Matrix')
         END IF
         write(iout,*) 'Opening File to Write Full Cholesky Decomposition'
         file_loc=File_Directory(4)(1:len_dir(4))//'/'//'Cholesky_Decomposition'
         Call pakstr(file_loc,len(2))
         Call IOsys('open cholesky_decomposition as new',0,0,0,file_loc(1:len(2)))
         Call IOsys('write real cholesky_factor to cholesky_decomposition',       &
                     2*tri_size,mat_z%upper,0,' ')
         write(iout,*) 'Closing File Holding Full Cholesky Decomposition'
      END IF
      count = 0
      DO i=1,n3d
         READ(50) rinput(1:i)
         ham(count+1:count+i) = rinput(1:i)
         count = count + i
      END DO
      Call scale_matrix(ham,smallest,largest,count,'hamiltonian')
      IF(print_cc) THEN          
         local_title='Input Hamiltonian Matrix'
         write(iout,1) local_title
         count = 0
         DO i=1,n3d
            write(iout,2) i
            write(iout,3) (ham(j), j=count+1, count + i )
            count = count + i
         END DO       
      END IF
  END IF
  DEALLOCATE(rinput)
1 Format(a80)
2 Format(/,5x,'Row = ',i4)
3 Format( (15x,6f10.5) )
! ***********************************************************************
  END SUBROUTINE Input_Full_Matrices_z
!***********************************************************************
!***********************************************************************
!deck Reformat_Matrices_d
!***begin prologue     Reformat_Matrices_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose
!***
!***references
!***routines called
!***end prologue       Reformat_Matrices_d
!
  SUBROUTINE Reformat_Matrices_d (ham,s)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:)               :: ham
  REAL(idp), DIMENSION(:)               :: s
  REAL(idp), DIMENSION(:), ALLOCATABLE  :: h_mat
  REAL(idp), DIMENSION(:), ALLOCATABLE  :: s_mat
  CHARACTER(LEN=3)                      :: itoc
  INTEGER                               :: i
  INTEGER                               :: j
!
  ALLOCATE( h_mat(new_tri_size) )
  Call New_Matrix(ham,h_mat)
  DEALLOCATE( h_mat )
  IF(print_cc) THEN          
     count = 0
     local_title='Reformatted Input Hamiltonian Matrix'
     write(iout,1) local_title
     DO i=1,new_size
        write(iout,2) i
        write(iout,3) (ham(j), j=count+1, count + i )
        count = count + i
     END DO       
  END IF
  IF (non_orth) THEN
      ALLOCATE( s_mat(new_tri_size) )
      Call New_Matrix(s,s_mat)
      DEALLOCATE( s_mat )
      IF(print_cc) THEN          
         count = 0
         local_title='Reformatted Input Overlap Matrix'
         write(iout,1) local_title
         DO i=1,new_size
            write(iout,2) i
            write(iout,3) (s(j), j=count+1, count + i )
            count = count + i
         END DO       
      END IF
  END IF
1 Format(a80)
2 Format(/,5x,'Row = ',i4)
3 Format( (15x,5f10.5) )
!***********************************************************************
  END SUBROUTINE Reformat_Matrices_d
!***********************************************************************
!***********************************************************************
!deck Reformat_Matrices_z
!***begin prologue     Reformat_Matrices_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose
!***
!***references
!***routines called
!***end prologue       Reformat_Matrices_z
!
  SUBROUTINE Reformat_Matrices_z (ham,s)
  IMPLICIT NONE
  COMPLEX(idp), DIMENSION(:)               :: ham
  COMPLEX(idp), DIMENSION(:)               :: s
  COMPLEX(idp), DIMENSION(:), ALLOCATABLE  :: h_mat
  COMPLEX(idp), DIMENSION(:), ALLOCATABLE  :: s_mat
  CHARACTER(LEN=3)                         :: itoc
  INTEGER                                  :: i
  INTEGER                                  :: j
!
  ALLOCATE( h_mat(new_tri_size) )
  Call New_Matrix (ham,h_mat)
  DEALLOCATE( h_mat)
  IF(print_cc) THEN          
     count = 0
     local_title='Reformatted Input Hamiltonian Matrix'
     write(iout,1) local_title
     DO i=1,new_size
        write(iout,2) i
        write(iout,3) (ham(j), j=count+1, count + i )
        count = count + i
     END DO       
  END IF
  IF (non_orth) THEN
      ALLOCATE( s_mat(new_tri_size) )
      Call New_Matrix (s,s_mat)
      DEALLOCATE( s_mat )
      IF(print_cc) THEN          
         count = 0
         local_title='Reformatted Input Overlap Matrix'
         write(iout,1) local_title
         DO i=1,new_size
            write(iout,2) i
            write(iout,3) (s(j), j=count+1, count + i )
            count = count + i
         END DO       
      END IF
  END IF
1 Format(a80)
2 Format(/,5x,'Row = ',i4)
3 Format( (15x,6f10.5) )
!***********************************************************************
  END SUBROUTINE Reformat_Matrices_z
!***********************************************************************
!***********************************************************************
!deck Input_Channel_Matrices_d
!***begin prologue     Input_Channel_Matrices_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Take the triangle of a real symmetric matrix and reformat it into
!***                   a set of channel submatrices.  The original matrix is deallocated.
!***references
!***routines called
!***end prologue       Input_Channel_Matrices_d
!
  SUBROUTINE Input_Channel_Matrices_d (ham,s)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:)               :: ham
  REAL(idp), DIMENSION(:)               :: s
  REAL(idp), DIMENSION(:), ALLOCATABLE  :: h_mat
  REAL(idp), DIMENSION(:), ALLOCATABLE  :: s_mat
  CHARACTER(LEN=3)                      :: itoc
  INTEGER                               :: i
  INTEGER                               :: j
!
  ALLOCATE( channel_mat(number_of_channels,number_of_channels) )
  DO ic=1, number_of_channels
     DO jc=1, number_of_channels
        ALLOCATE(channel_mat(ic,jc)%mat_d%channel_h_matrix(channel_size,channel_size))
     END DO
  END DO
  ALLOCATE( h_mat(new_tri_size) )
  Call Channel_Sub_Matrix(ham,h_mat,'hamiltonian')
  IF(print_cc) THEN
     DO ic = 1, number_of_channels
         DO jc = 1, ic
            Call Print_Matrix(type_real_matrix,channel_mat(ic,jc)%mat_d%channel_h_matrix,    &
                             new_size, new_size, title='Hamiltonian Matrix for channel i = ' &
                            //itoc(ic)//' channel j = '//itoc(jc) )
        END DO
     END DO
  END IF
  DEALLOCATE( h_mat )
  IF (non_orth) THEN
      DO ic=1, number_of_channels
         DO jc=1, number_of_channels
            ALLOCATE(                                                                        &
            channel_mat(ic,jc)%mat_d%channel_s_matrix(channel_size,channel_size))
         END DO
      END DO
      ALLOCATE( s_mat(new_tri_size) )
      Call Channel_Sub_Matrix(s,s_mat,'overlap')
      DEALLOCATE( s_mat )
      IF(print_cc) THEN
          DO ic = 1, number_of_channels
             DO jc = 1, ic
                Call Print_Matrix(type_real_matrix,channel_mat(ic,jc)%mat_d%channel_s_matrix, &
                                 new_size, new_size, title='Overlap Matrix for channel i = '  &
                                 //itoc(ic)//' channel j = '//itoc(jc) )
             END DO
          END DO
      END IF
  END IF
!***********************************************************************
  END SUBROUTINE Input_Channel_Matrices_d
!***********************************************************************
!***********************************************************************
!deck Input_Channel_Matrices_z
!***begin prologue     Input_Channel_Matrices_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Take the triangle of a Hermitian matrix and reformat it into
!***                   a set of channel submatrices.  The original matrix is deallocated.
!***
!***references
!***routines called
!***end prologue       Input_Channel_Matrices_z
!
  SUBROUTINE Input_Channel_Matrices_z (ham,s)
  IMPLICIT NONE
  COMPLEX(idp), DIMENSION(:)               :: ham
  COMPLEX(idp), DIMENSION(:)               :: s
  COMPLEX(idp), DIMENSION(:), ALLOCATABLE  :: h_mat
  COMPLEX(idp), DIMENSION(:), ALLOCATABLE  :: s_mat
  CHARACTER(LEN=3)                         :: itoc
  INTEGER                                  :: i
  INTEGER                                  :: j
!
  ALLOCATE( channel_mat(number_of_channels,number_of_channels) )
  DO ic=1, number_of_channels
     DO jc=1, number_of_channels
        ALLOCATE(channel_mat(ic,jc)%mat_z%channel_h_matrix(channel_size,channel_size))
     END DO
  END DO
  ALLOCATE( h_mat(new_tri_size) )
  Call Channel_Sub_Matrix(ham,h_mat,'hamiltonian')
  DEALLOCATE( h_mat)
  IF(print_channel_matrices) THEN
     DO ic = 1, number_of_channels
         DO jc = 1, ic
            Call Print_Matrix(type_complex_matrix,channel_mat(ic,jc)%mat_z%channel_h_matrix,   &
                             new_size, new_size, title='Hamiltonian Matrix for channel i = '   &
                            //itoc(ic)//' channel j = '//itoc(jc) )
        END DO
     END DO
  END IF
  IF (non_orth) THEN
      DO ic=1, number_of_channels
         DO jc=1, number_of_channels
            ALLOCATE(channel_mat(ic,jc)%mat_z%channel_s_matrix(channel_size,channel_size))
         END DO
      END DO
      ALLOCATE( s_mat(new_tri_size) )
      Call Channel_Sub_Matrix(s,s_mat,'overlap')
      DEALLOCATE( s_mat )
      IF(print_channel_matrices) THEN
          DO ic = 1, number_of_channels
             DO jc = 1, ic
                Call Print_Matrix(type_complex_matrix,channel_mat(ic,jc)%mat_z%channel_s_matrix,      &
                                  new_size, new_size, title='Overlap Matrix for channel i = '         &
                               //itoc(ic)//' channel j = '//itoc(jc) )
             END DO
          END DO
      END IF
  END IF
!***********************************************************************
  END SUBROUTINE Input_Channel_Matrices_z
!***********************************************************************
!***********************************************************************
!deck Channel_Sub_Matrix_d
!***begin prologue     Channel_Sub_Matrix_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose
!***
!***references
!***routines called
!***end prologue       Channel_Sub_Matrix_d
!
  SUBROUTINE Channel_Sub_Matrix_d (old_matrix,new_matrix,matrix_type)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:)               :: old_matrix
  REAL(idp), DIMENSION(:)               :: new_matrix
  CHARACTER(LEN=*)                      :: matrix_type
  CHARACTER(LEN=3)                      :: itoc
  INTEGER                               :: i
  INTEGER                               :: j
  IF (matrix_type =='hamiltonian') THEN
!
     new_matrix_count_i = 0
     DO ic=1, number_of_channels
        new_matrix_count_j = 0
        DO jc = 1, ic
           IF ( ic == jc) THEN
                number_of_splines_i = 0
                ii = new_matrix_count_i
                DO is = first_spline, last_spline
                   number_of_splines_i = number_of_splines_i + 1
                   i = channel_buf(ic,is)
                   ii = ii + 1
                   new_ind = ii * ( ii - 1 ) / 2 
                   old_ind =  i * ( i-1 ) / 2
                   number_of_splines_j = 0
                   jj = new_matrix_count_j
                   DO js = first_spline, is
                      number_of_splines_j = number_of_splines_j + 1
                      jj = jj + 1
                      j = channel_buf(jc,js)
                      old_ij = old_ind + j
                      new_ij = new_ind + jj
                      channel_mat(ic,jc)%mat_d%channel_h_matrix(number_of_splines_i, &
                                                                number_of_splines_j) &
                                            = old_matrix( old_ij )
                      new_matrix( new_ij )  = old_matrix( old_ij )
                   END DO
                END DO
           ELSE 
                number_of_splines_i = 0
                ii = new_matrix_count_i
                DO is = first_spline, last_spline
                   number_of_splines_i = number_of_splines_i + 1
                   i = channel_buf(ic,is)
                   ii = ii + 1
                   new_ind = ii * ( ii -1 ) / 2 
                   old_ind =  i * ( i-1 ) / 2
                   number_of_splines_j = 0
                   jj = new_matrix_count_j
                   DO js = first_spline, last_spline
                      number_of_splines_j = number_of_splines_j + 1
                      jj = jj + 1
                      j = channel_buf(jc,js)
                      old_ij = old_ind + j
                      new_ij = new_ind + jj
                      channel_mat(ic,jc)%mat_d%channel_h_matrix(number_of_splines_i, &
                                                                number_of_splines_j) &
                                            = old_matrix( old_ij )
                      new_matrix( new_ij )  = old_matrix( old_ij )
                   END DO
                END DO
           END IF
           new_matrix_count_j = new_matrix_count_j + channel_size
        END DO
        new_matrix_count_i = new_matrix_count_i + channel_size
     END DO
  ELSE IF (matrix_type == 'overlap') THEN
!
     new_matrix_count_i = 0
     DO ic=1, number_of_channels
        new_matrix_count_j = 0
        DO jc = 1, ic
           IF ( ic == jc) THEN
                number_of_splines_i = 0
                ii = new_matrix_count_i
                DO is = first_spline, last_spline
                   number_of_splines_i = number_of_splines_i + 1
                   i = channel_buf(ic,is)
                   ii = ii + 1
                   new_ind = ii * ( ii - 1 ) / 2 
                   old_ind =  i * ( i-1 ) / 2
                   number_of_splines_j = 0
                   jj = new_matrix_count_j
                   DO js = first_spline, is
                      number_of_splines_j = number_of_splines_j + 1
                      jj = jj + 1
                      j = channel_buf(jc,js)
                      old_ij = old_ind + j
                      new_ij = new_ind + jj
                      channel_mat(ic,jc)%mat_d%channel_s_matrix(number_of_splines_i, &
                                                                number_of_splines_j) &
                                            = old_matrix( old_ij )
                      new_matrix( new_ij )  = old_matrix( old_ij )
                   END DO
                END DO
           ELSE 
                number_of_splines_i = 0
                ii = new_matrix_count_i
                DO is = first_spline, last_spline
                   number_of_splines_i = number_of_splines_i + 1
                   i = channel_buf(ic,is)
                   ii = ii + 1
                   new_ind = ii * ( ii -1 ) / 2 
                   old_ind =  i * ( i-1 ) / 2
                   number_of_splines_j = 0
                   jj = new_matrix_count_j
                   DO js = first_spline, last_spline
                      number_of_splines_j = number_of_splines_j + 1
                      jj = jj + 1
                      j = channel_buf(jc,js)
                      old_ij = old_ind + j
                      new_ij = new_ind + jj
                      channel_mat(ic,jc)%mat_d%channel_s_matrix(number_of_splines_i, &
                                                                number_of_splines_j) &
                                            = old_matrix( old_ij )
                      new_matrix( new_ij )  = old_matrix( old_ij )
                   END DO
                END DO
           END IF
           new_matrix_count_j = new_matrix_count_j + channel_size
        END DO
        new_matrix_count_i = new_matrix_count_i + channel_size
     END DO
  END IF
  old_matrix(1:new_tri_size) = new_matrix(1:new_tri_size)
!***********************************************************************
  END SUBROUTINE Channel_Sub_Matrix_d
!***********************************************************************
!deck Channel_Sub_Matrix_z
!***begin prologue     Channel_Sub_Matrix_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose
!***
!***references
!***routines called
!***end prologue       Channel_Sub_Matrix_z
!
  SUBROUTINE Channel_Sub_Matrix_z (old_matrix,new_matrix,matrix_type)
  IMPLICIT NONE
  COMPLEX(idp), DIMENSION(:)            :: old_matrix
  COMPLEX(idp), DIMENSION(:)            :: new_matrix
  CHARACTER(LEN=*)                      :: matrix_type
  CHARACTER(LEN=3)                      :: itoc
  INTEGER                               :: i
  INTEGER                               :: j
  IF (matrix_type =='hamiltonian') THEN
!
     new_matrix_count_i = 0
     DO ic=1, number_of_channels
        new_matrix_count_j = 0
        DO jc = 1, ic
           IF ( ic == jc) THEN
                number_of_splines_i = 0
                ii = new_matrix_count_i
                DO is = first_spline, last_spline
                   number_of_splines_i = number_of_splines_i + 1
                   i = channel_buf(ic,is)
                   ii = ii + 1
                   new_ind = ii * ( ii - 1 ) / 2 
                   old_ind =  i * ( i-1 ) / 2
                   number_of_splines_j = 0
                   jj = new_matrix_count_j
                   DO js = first_spline, is
                      number_of_splines_j = number_of_splines_j + 1
                      jj = jj + 1
                      j = channel_buf(jc,js)
                      old_ij = old_ind + j
                      new_ij = new_ind + jj
                      channel_mat(ic,jc)%mat_z%channel_h_matrix(number_of_splines_i, &
                                                                number_of_splines_j) &
                                            = old_matrix( old_ij )
                      new_matrix( new_ij )  = old_matrix( old_ij )
                   END DO
                END DO
           ELSE 
                number_of_splines_i = 0
                ii = new_matrix_count_i
                DO is = first_spline, last_spline
                   number_of_splines_i = number_of_splines_i + 1
                   i = channel_buf(ic,is)
                   ii = ii + 1
                   new_ind = ii * ( ii -1 ) / 2 
                   old_ind =  i * ( i-1 ) / 2
                   number_of_splines_j = 0
                   jj = new_matrix_count_j
                   DO js = first_spline, last_spline
                      number_of_splines_j = number_of_splines_j + 1
                      jj = jj + 1
                      j = channel_buf(jc,js)
                      old_ij = old_ind + j
                      new_ij = new_ind + jj
                      channel_mat(ic,jc)%mat_z%channel_h_matrix(number_of_splines_i, &
                                                                number_of_splines_j) &
                                            = old_matrix( old_ij )
                      new_matrix( new_ij )  = old_matrix( old_ij )
                   END DO
                END DO
           END IF
           new_matrix_count_j = new_matrix_count_j + channel_size
        END DO
        new_matrix_count_i = new_matrix_count_i + channel_size
     END DO
  ELSE IF (matrix_type == 'overlap') THEN
!
     new_matrix_count_i = 0
     DO ic=1, number_of_channels
        new_matrix_count_j = 0
        DO jc = 1, ic
           IF ( ic == jc) THEN
                number_of_splines_i = 0
                ii = new_matrix_count_i
                DO is = first_spline, last_spline
                   number_of_splines_i = number_of_splines_i + 1
                   i = channel_buf(ic,is)
                   ii = ii + 1
                   new_ind = ii * ( ii - 1 ) / 2 
                   old_ind =  i * ( i-1 ) / 2
                   number_of_splines_j = 0
                   jj = new_matrix_count_j
                   DO js = first_spline, is
                      number_of_splines_j = number_of_splines_j + 1
                      jj = jj + 1
                      j = channel_buf(jc,js)
                      old_ij = old_ind + j
                      new_ij = new_ind + jj
                      channel_mat(ic,jc)%mat_z%channel_s_matrix(number_of_splines_i, &
                                                                number_of_splines_j) &
                                            = old_matrix( old_ij )
                      new_matrix( new_ij )  = old_matrix( old_ij )
                   END DO
                END DO
           ELSE 
                number_of_splines_i = 0
                ii = new_matrix_count_i
                DO is = first_spline, last_spline
                   number_of_splines_i = number_of_splines_i + 1
                   i = channel_buf(ic,is)
                   ii = ii + 1
                   new_ind = ii * ( ii -1 ) / 2 
                   old_ind =  i * ( i-1 ) / 2
                   number_of_splines_j = 0
                   jj = new_matrix_count_j
                   DO js = first_spline, last_spline
                      number_of_splines_j = number_of_splines_j + 1
                      jj = jj + 1
                      j = channel_buf(jc,js)
                      old_ij = old_ind + j
                      new_ij = new_ind + jj
                      channel_mat(ic,jc)%mat_z%channel_s_matrix(number_of_splines_i, &
                                                            number_of_splines_j) &
                                            = old_matrix( old_ij )
                      new_matrix( new_ij )  = old_matrix( old_ij )
                   END DO
                END DO
           END IF
           new_matrix_count_j = new_matrix_count_j + channel_size
        END DO
        new_matrix_count_i = new_matrix_count_i + channel_size
     END DO
  END IF
  old_matrix(1:new_tri_size) = new_matrix(1:new_tri_size)
!***********************************************************************
  END SUBROUTINE Channel_Sub_Matrix_z
!***********************************************************************
!***********************************************************************
!Compute_Eigenstates_d
!***begin prologue     Compute_Eigenstates_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Compute_Eigenstates_d
!
  SUBROUTINE Compute_Eigenstates_d (ham,s)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:)                 :: ham
  REAL(idp), DIMENSION(:)                 :: s
  REAL(idp), DIMENSION(:,:), ALLOCATABLE  :: h_mat
  REAL(idp), DIMENSION(:,:), ALLOCATABLE  :: s_mat
  REAL(idp), DIMENSION(:),   ALLOCATABLE  :: eig
  REAL(idp), DIMENSION(:),   ALLOCATABLE  :: work
  CHARACTER (LEN=8)                       :: itoc
  CHARACTER (LEN=8)                       :: kewrd
  INTEGER                                 :: i
  INTEGER                                 :: j
!
   ALLOCATE( h_mat(1:n3d,1:n3d), eig(1:n3d), work(1:lwork) )
   count = 0
   DO i = 1, n3d
      DO j = 1, i
         count = count + 1
         h_mat(i,j) = ham(count) 
         h_mat(j,i) = ham(count) 
      END DO
   END DO
   IF(non_orth) THEN      
      ALLOCATE( s_mat(1:n3d,1:n3d) )
      count = 0
      DO i = 1,n3d
         DO j=1,i
            count = count + 1
            s_mat(i,j) = s(count) 
            s_mat(j,i) = s(count) 
         END DO
      END DO
   END IF
   IF(non_orth) THEN
      call dsygv(1,'v','u',n3d,h_mat,n3d,s_mat,n3d,eig,work,lwork,info)
   ELSE
      call dsyev('v','u',n3d,h_mat,n3d,eig,work,lwork,info)
   END IF
   local_title='Exact Eigenvalues'
   call Print_Matrix(type_real_vector,eig,title='Exact Eigenvalues')
   IF(print_cc) THEN
      call Print_Matrix(type_real_matrix,h_mat,n3d,n3d,title='Exact Eigenvectors')
   END IF
   IF (input_output_vector == 'output') THEN
       kewrd='real'
       write(iout,*) 'Writing '//kewrd//' data to eigenstates'
       write(20) kewrd
       DO i=1,n3d 
          write(20) (h_mat(j,i), j=1,n3d)
       END DO
       CLOSE(20)
   END IF
   DEALLOCATE( h_mat, eig, work )
   IF(non_orth) THEN      
      DEALLOCATE( s_mat )
   END IF
!***********************************************************************
  END SUBROUTINE Compute_Eigenstates_d
!***********************************************************************
!***********************************************************************
!deck Compute_Eigenstates_z
!***begin prologue     Compute_Eigenstates_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Compute_Eigenstates_z
!
  SUBROUTINE Compute_Eigenstates_z(ham,s)
  IMPLICIT NONE
  COMPLEX(idp), DIMENSION(:)                 :: ham 
  COMPLEX(idp), DIMENSION(:)                 :: s 
  COMPLEX(idp), DIMENSION(:,:), ALLOCATABLE  :: h_mat
  COMPLEX(idp), DIMENSION(:,:), ALLOCATABLE  :: s_mat
  REAL(idp),    DIMENSION(:),   ALLOCATABLE  :: eig
  COMPLEX(idp), DIMENSION(:),   ALLOCATABLE  :: work
  REAL(idp),    DIMENSION(:),   ALLOCATABLE  :: rwork
  CHARACTER (LEN=8)                          :: itoc
  CHARACTER (LEN=8)                          :: kewrd
  INTEGER                                    :: i
  INTEGER                                    :: j
   ALLOCATE( h_mat(1:n3d,1:n3d), eig(1:n3d), work(1:lwork), rwork(1:lwork) )
   count = 0
   DO i = 1, n3d
      DO j = 1, i
         count = count + 1
         h_mat(i,j) = ham(count)  
         h_mat(j,i) = conjg(ham(count))
      END DO
   END DO  
   IF(non_orth) THEN      
      ALLOCATE( s_mat(1:n3d,1:n3d) )
      count = 0
      DO i = 1, n3d
         DO j = 1, i
            count = count + 1
            s_mat(i,j) = s(count)  
            s_mat(j,i) = conjg(s(count))
         END DO
      END DO  
   END IF
  IF(non_orth) THEN
     call zhegv(1,'v','u',n3d,h_mat,n3d,s_mat,n3d,eig,work,lwork,rwork,info)
  ELSE
     call zheev('v','u',n3d,h_mat,n3d,eig,work,lwork,rwork,info)
  END IF
  local_title='Exact Eigenvalues'
  call Print_Matrix(type_real_vector,eig,title='Exact Eigenvalues')
  IF(print_cc) THEN
     local_title='Exact Eigenvectors'
     Call Print_Matrix(type_complex_matrix,h_mat,n3d,n3d,title='Exact Eigenvectors')
  END IF
  IF (input_output_vector == 'output') THEN
      kewrd='complex'
      write(iout,*) 'Writing '//kewrd//' data to eigenstates'
      write(20) kewrd
      DO i=1,n3d
         write(20) (h_mat(j,i), j=1,n3d)
      END DO
      CLOSE(20)
  END IF
  DEALLOCATE( h_mat, eig, work, rwork )
  IF(non_orth) THEN      
     DEALLOCATE( s_mat)
  END IF
! ***********************************************************************
  END SUBROUTINE Compute_Eigenstates_z
!***********************************************************************
!***********************************************************************
!deck New_Matrix_d
!***begin prologue     New_Matrix_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose
!***
!***references
!***routines called
!***end prologue       New_Matrix_d
!
  SUBROUTINE New_Matrix_d (old_matrix,new_matrix)
  IMPLICIT NONE
  REAL(idp), DIMENSION(:)               :: old_matrix
  REAL(idp), DIMENSION(:)               :: new_matrix
  CHARACTER(LEN=3)                      :: itoc
  INTEGER                               :: i
  INTEGER                               :: j
  new_matrix_count_i = 0
  DO ic=1, number_of_channels
     new_matrix_count_j = 0
     DO jc = 1, ic
        IF ( ic == jc) THEN
             number_of_splines_i = 0
             ii = new_matrix_count_i
             DO is = first_spline, last_spline
                number_of_splines_i = number_of_splines_i + 1
                i = channel_buf(ic,is)
                ii = ii + 1
                new_ind = ii * ( ii - 1 ) / 2 
                old_ind =  i * ( i-1 ) / 2
                number_of_splines_j = 0
                jj = new_matrix_count_j
                DO js = first_spline, is
                   number_of_splines_j = number_of_splines_j + 1
                   jj = jj + 1
                   j = channel_buf(jc,js)
                   old_ij = old_ind + j
                   new_ij = new_ind + jj
                   new_matrix( new_ij )  = old_matrix( old_ij )
                END DO
             END DO
        ELSE 
             number_of_splines_i = 0
             ii = new_matrix_count_i
             DO is = first_spline, last_spline
                number_of_splines_i = number_of_splines_i + 1
                i = channel_buf(ic,is)
                ii = ii + 1
                new_ind = ii * ( ii -1 ) / 2 
                old_ind =  i * ( i-1 ) / 2
                number_of_splines_j = 0
                jj = new_matrix_count_j
                DO js = first_spline, last_spline
                   number_of_splines_j = number_of_splines_j + 1
                   jj = jj + 1
                   j = channel_buf(jc,js)
                   old_ij = old_ind + j
                   new_ij = new_ind + jj
                   new_matrix( new_ij )  = old_matrix( old_ij )
                END DO
             END DO
        END IF
        new_matrix_count_j = new_matrix_count_j + channel_size
     END DO
     new_matrix_count_i = new_matrix_count_i + channel_size
  END DO
  old_matrix(1:new_tri_size) = new_matrix(1:new_tri_size)
!***********************************************************************
  END SUBROUTINE New_Matrix_d
!***********************************************************************
!***********************************************************************
!deck New_Matrix_z
!***begin prologue     New_Matrix_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose
!***
!***references
!***routines called
!***end prologue       New_Matrix_z
!
  SUBROUTINE New_Matrix_z (old_matrix,new_matrix)
  IMPLICIT NONE
  COMPLEX(idp), DIMENSION(:)            :: old_matrix
  COMPLEX(idp), DIMENSION(:)            :: new_matrix
  CHARACTER(LEN=3)                      :: itoc
  INTEGER                               :: i
  INTEGER                               :: j
  new_matrix_count_i = 0
  DO ic=1, number_of_channels
     new_matrix_count_j = 0
     DO jc = 1, ic
        IF ( ic == jc) THEN
             number_of_splines_i = 0
             ii = new_matrix_count_i
             DO is = first_spline, last_spline
                number_of_splines_i = number_of_splines_i + 1
                i = channel_buf(ic,is)
                ii = ii + 1
                new_ind = ii * ( ii - 1 ) / 2 
                old_ind =  i * ( i-1 ) / 2
                number_of_splines_j = 0
                jj = new_matrix_count_j
                DO js = first_spline, is
                   number_of_splines_j = number_of_splines_j + 1
                   jj = jj + 1
                   j = channel_buf(jc,js)
                   old_ij = old_ind + j
                   new_ij = new_ind + jj
                   new_matrix( new_ij )  = old_matrix( old_ij )
                END DO
             END DO
        ELSE 
             number_of_splines_i = 0
             ii = new_matrix_count_i
             DO is = first_spline, last_spline
                number_of_splines_i = number_of_splines_i + 1
                i = channel_buf(ic,is)
                ii = ii + 1
                new_ind = ii * ( ii -1 ) / 2 
                old_ind =  i * ( i-1 ) / 2
                number_of_splines_j = 0
                jj = new_matrix_count_j
                DO js = first_spline, last_spline
                   number_of_splines_j = number_of_splines_j + 1
                   jj = jj + 1
                   j = channel_buf(jc,js)
                   old_ij = old_ind + j
                   new_ij = new_ind + jj
                   new_matrix( new_ij )  = old_matrix( old_ij )
                END DO
             END DO
        END IF
        new_matrix_count_j = new_matrix_count_j + channel_size
     END DO
     new_matrix_count_i = new_matrix_count_i + channel_size
  END DO
  old_matrix(1:new_tri_size) = new_matrix(1:new_tri_size)
!***********************************************************************
  END SUBROUTINE New_Matrix_z
!***********************************************************************
!***********************************************************************
!Compute_Inverse_d
!***begin prologue     Compute_Inverse_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Compute_Inverse_d
!
  SUBROUTINE Compute_Inverse_d(s)
  USE full_matrix_vector_iteration_module
  IMPLICIT NONE
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  CHARACTER(LEN=1)                      :: rowlab
  CHARACTER(LEN=1)                      :: collab
  LOGICAL                               :: dollar
  LOGICAL                               :: logkey
  INTEGER                               :: intkey
  INTEGER                               :: tri
  INTEGER                               :: number_of_right_hand_sides
  INTEGER                               :: nrhs
  INTEGER                               :: size
  REAL(idp)                             :: fpkey
  REAL(idp)                             :: omega
  REAL(idp), DIMENSION(:)               :: s
  REAL(idp), DIMENSION(:), ALLOCATABLE  :: s_vec_in
  REAL(idp), DIMENSION(:), ALLOCATABLE  :: s_vec_out
  REAL(idp), DIMENSION(:), ALLOCATABLE  :: rhs
!
  IF ( dollar('$iteration_parameters',data_card,pass_data,inp) ) then
       iteration_method = chrkey(data_card,'iteration_method','jacobi',' ')
       size = intkey(data_card,'size_of_matrix',2,' ')
       number_of_right_hand_sides = intkey(data_card,'number_of_right_hand_sides',  &
                                           1,' ')
       n_iter = intkey(data_card,'number_of_iterations',size,' ')
       tol = fpkey(data_card,'iteration_tolerance',1.d-10,' ')
       omega = fpkey(data_card,'relaxation_parameter',1.d0,' ')
       iterative_print = logkey(data_card,'iterative_print',.false.,' ')
       tri = size * ( size + 1 ) / 2
!       ALLOCATE( s_matrix(1:tri) )
!       Call fparr(data_card,'s_matrix',s_matrix,tri,' ')
  END IF
  size = n3d
  ALLOCATE( s_vec_in(1:size), s_vec_out(1:size), rhs(1:size) )
  Write(iout,1) iteration_method, size, n_iter, tol
  DO nrhs = 1, number_of_right_hand_sides
     write(iout,*) 'right hand side = ',nrhs
     call random_number(rhs(1:size))
!     rhs_d(1) = 1.d0
!     rhs_d(10) = 1.d0
!     rhs_d(20) = 1.d0
!     rhs_d(100) = 1.d0
     IF (iteration_method == 'jacobi') THEN
         Call jacobi_iteration(s,s_vec_in,s_vec_out,rhs)  
     ELSE IF (iteration_method == 'gauss_seidel') THEN
         Call gauss_seidel_iteration(s,s_vec_in,s_vec_out,rhs,omega)  
     END IF
  END DO
1 Format(/,5x,'iteration method     = ',a16,2x,'Matrix Size = ',i5,         &
         /,5x,'number of Iterations = ',i5,2x,'Tolerance    = ',e15.8)
2 Format(a80)
3 Format(/,5x,'Row = ',i4)
4 Format( (15x,5e15.8) )
  DEALLOCATE( s_vec_in, s_vec_out, rhs )
!***********************************************************************
  END SUBROUTINE Compute_Inverse_d
!***********************************************************************
!***********************************************************************
!Compute_Inverse_z
!***begin prologue     Compute_Inverse_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Compute_Inverse_z
!
  SUBROUTINE Compute_Inverse_z (s)
  USE full_matrix_vector_iteration_module
  IMPLICIT NONE
  CHARACTER(LEN=80)                        :: chrkey
  CHARACTER (LEN=8)                        :: itoc
  LOGICAL                                  :: dollar
  LOGICAL                                  :: logkey
  INTEGER                                  :: intkey
  INTEGER                                  :: size
  REAL(idp)                                :: fpkey
  COMPLEX(idp), DIMENSION(:)               :: s
  COMPLEX(idp), DIMENSION(:), ALLOCATABLE  :: s_vec_in
  COMPLEX(idp), DIMENSION(:), ALLOCATABLE  :: s_vec_out
  COMPLEX(idp), DIMENSION(:), ALLOCATABLE  :: rhs
  REAL(idp)                                :: t_ran
  REAL(idp)                                :: u_ran
  REAL(idp)                                :: omega
  IF ( dollar('$iteration_parameters',data_card,pass_data,inp) ) then
       n_iter = intkey(data_card,'number_of_iterations',n3d,' ')
       tol = fpkey(data_card,'iteration_tonlerance',1.d-10,' ')
  END IF
  ALLOCATE( s_vec_in(1:n3d), s_vec_out(1:n3d), rhs(1:n3d) )
  size = n3d
  DO i = 1, n3d
     call random_number(t_ran)
     call random_number(u_ran)
     rhs(i) = cmplx(t_ran,u_ran)
   END DO
  IF (iteration_method == 'jacobi') THEN
      Call jacobi_iteration(s,s_vec_in,s_vec_out,rhs)  
  ELSE IF (iteration_method == 'gauss_seidel') THEN
      Call gauss_seidel_iteration(s,s_vec_in,s_vec_out,rhs,omega)  
  END IF
   DEALLOCATE( s_vec_in, s_vec_out, rhs )
!***********************************************************************
  END SUBROUTINE Compute_Inverse_z
!***********************************************************************
!***********************************************************************
!***begin prologue     scale_matrix__d
!***date written       000710   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           one-dim
!***author             schneider, barry (nsf)
!***source             
!***purpose            Pack non zero hamiltonian matrix elements and indices. 
!***                   What is stored are the non zero elements of each column
!***                   of the upper triangle.  The array packed_columns has
!***                   these.  The row index is stored is row_index and the
!***                   number of non zero elements in the column in non_zero_column.
!***                   Storage is a11, a12, a22, a13, a23, a33 etc.
!***references         
!
!***routines called    
!***end prologue       scale_matrix_d
  Subroutine scale_matrix_d(matrix,smallest,largest,n,type)
  IMPLICIT NONE
  REAL(idp),   DIMENSION(:)    :: matrix
  REAL(idp)                    :: smallest
  REAL(idp)                    :: largest
  INTEGER                      :: n
  INTEGER                      :: i
  INTEGER                      :: j
  INTEGER                      :: num
  INTEGER                      :: number
  INTEGER                      :: count
  LOGICAL                      :: prn
  CHARACTER(LEN=*)             :: type
  CHARACTER(LEN=80)            :: title
!
  smallest = abs(matrix(1))
  largest  = smallest
  DO i=1,n
     smallest = min(smallest,abs(matrix(i)))
     largest  = max(largest,abs(matrix(i)))
  END DO
!
  WRITE(iout,1) type
  WRITE(iout,2) smallest, largest
1 FORMAT(/,10x,'Matrix Type = ', a16)
2 FORMAT(/,10x,'Absolute Smallest Element =', e18.8,5x, 'Absolute Largest Element =',e15.8)
!
END SUBROUTINE scale_matrix_d
!***********************************************************************
!***********************************************************************
!***begin prologue     scale_matrix_z
!***date written       000710   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           one-dim
!***author             schneider, barry (nsf)
!***source             
!***purpose            Pack non zero hamiltonian matrix elements and indices. 
!***                   What is stored are the non zero elements of each column
!***                   of the upper triangle.  The array packed_columns has
!***                   these.  The row index is stored is row_index and the
!***                   number of non zero elements in the column in non_zero_column.
!***                   Storage is a11, a12, a22, a13, a23, a33 etc.
!***references         
!
!***routines called    
!***end prologue       scale_matrix_z
  Subroutine scale_matrix_z(matrix,smallest,largest,n,type)
  IMPLICIT NONE
  COMPLEX(idp),   DIMENSION(:) :: matrix
  REAL(idp)                    :: smallest
  REAL(idp)                    :: largest
  INTEGER                      :: n
  INTEGER                      :: i
  INTEGER                      :: j
  INTEGER                      :: num
  INTEGER                      :: number
  INTEGER                      :: count
  LOGICAL                      :: prn
  CHARACTER(LEN=*)             :: type
  CHARACTER(LEN=80)            :: title
!
  smallest = abs(matrix(1))
  largest  = smallest
  DO i=1,n
     smallest = min(smallest,abs(matrix(i)))
     largest  = max(largest,abs(matrix(i)))
  END DO
!
  WRITE(iout,1) type
  WRITE(iout,2) smallest, largest
1 FORMAT(/,10x,'Matrix Type = ', a16)
2 FORMAT(/,10x,'Absolute Smallest Element =', e18.8,5x, 'Absolute Largest Element =',e15.8)
!
END SUBROUTINE scale_matrix_z
!***********************************************************************
!***********************************************************************
  END  MODULE Non_Packed_Matrix_Module
!***********************************************************************
!***********************************************************************
