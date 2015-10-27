!***********************************************************************
                           MODULE B_Spline_Matrix_Module
!
                           USE Global_Time_Propagation_Module
                           USE Pack_Hamiltonian_Module
                           USE Iterative_Global,  ONLY : overlap_d, hamiltonian_d
                           USE Prop_Global,  ONLY : n3d
!
                           IMPLICIT NONE
            CHARACTER(LEN=1600)                  :: data_card
            CHARACTER(LEN=80)                    :: pass_data
            CHARACTER(LEN=80)                    :: local_title
!***********************************************************************
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Deck Read_Channel_B_Matrices
!***begin prologue     Read_Channel_B_Spline_Matrices
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Read in B spline matrices in Oleg's format
!***
!***references
!***routines called
!***end prologue       Read_Channel_B_Spline_Matrices
!
  SUBROUTINE Read_Channel_B_Spline_Matrices(file_directory, matrix_file )
  IMPLICIT NONE
  CHARACTER(LEN=*)                      :: file_directory
  CHARACTER(LEN=*)                      :: matrix_file
  INTEGER                               :: lenth
  INTEGER                               :: len
  INTEGER                               :: len_1
  INTEGER                               :: IOSTAT
!
  Call pakstr(file_directory,len)
  Call pakstr(matrix_file,len_1) 
  OPEN(UNIT=50,FILE=file_directory(1:len)//'/'//matrix_file(1:len_1),                 &
       ACCESS='sequential',FORM='unformatted',IOSTAT=IOSTAT,STATUS='old')
  Read(50) 
  Read(50)
  read(50)
  ALLOCATE(overlap_d(n3d,n3d))
  overlap_d = 0.d0
  Call Read_Mat(overlap_d,spline_array_j)
  ALLOCATE(hamiltonian_d(n3d,n3d))
  hamiltonian_d = 0.d0
  Call Read_Mat(hamiltonian_d,spline_array_j)
  stop
!***********************************************************************
  END SUBROUTINE Read_Channel_B_Spline_Matrices
!***********************************************************************
!***********************************************************************
!Deck Read_Mat
!***begin prologue     Read_Mat
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Read in B spline matrices in Oleg's format
!***
!***references
!***routines called
!***end prologue       Read_Mat
!
  SUBROUTINE Read_Mat(matrix,spline_array)
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:)                :: matrix
  INTEGER, DIMENSION(0:*)               :: spline_array
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: i_1
  INTEGER                               :: i_2
  INTEGER                               :: j_1
  INTEGER                               :: j_2
  INTEGER                               :: i_ch
  INTEGER                               :: j_ch
  INTEGER                               :: size_i
  INTEGER                               :: size_j
  CHARACTER(LEN=3)                      :: itoc
!
  DO
    Read(50) i_ch, j_ch
    IF(i_ch <= 0) THEN
       EXIT
    END IF
    Write(iout,*) 'Channel Block 1 = ', i_ch, 'Channel Block 2 = ', j_ch
    i_1 = spline_array(i_ch - 1 ) + 1
    i_2 = spline_array(i_ch)  
    j_1 = spline_array(j_ch - 1 ) + 1
    j_2 = spline_array(j_ch)  
    size_i = i_2 - i_1 + 1
    size_j = j_2 - j_1 + 1
    DO j =j_1, j_2
       Read(50) (matrix(i,j), i=i_1,i_2)
    END DO
    IF ( i_ch /= j_ch ) THEN
         DO j = j_1, j_2
            matrix(j,i_1:i_2) = matrix(i_1:i_2,j)
         END DO
    END IF
    IF(print_cc) THEN          
       local_title='Matrix Channel 1 = '//itoc(i_ch)//'Channel 2 = '//itoc(j_ch)
       call prntfmn(local_title,matrix(i_1,j_1),size_i,size_j,n3d,n3d,iout,'e')
    END IF
  END DO
!***********************************************************************
  END SUBROUTINE Read_Mat
!***********************************************************************
!***********************************************************************
  END  MODULE B_Spline_Matrix_Module
!***********************************************************************
!***********************************************************************
