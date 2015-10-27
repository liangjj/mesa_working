!***********************************************************************
                           MODULE Matrix_Module
                           USE dvrprop_global
                           USE Preconditioner_Module
                           USE Iterative_Global
                           USE Pack_Matrix_Module
!
                           IMPLICIT NONE


  CHARACTER (LEN=80)                       :: matrix_file
  CHARACTER (LEN=80)                       ::  matrix_type
  INTEGER                                  :: n_splines
  INTEGER                                  :: n_chan
  INTEGER                                  :: n_bound
  LOGICAL                                  :: ifeig=.false.
  LOGICAL                                  :: solve_only
  LOGICAL                                  :: read_only
!***********************************************************************
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          INTERFACE Compute_Solutions
                   MODULE PROCEDURE Compute_Solutions_d,                   &
                                    Compute_Solutions_Tri_d,               &
                                    Compute_Solutions_z,                   &
                                    Compute_Solutions_Tri_z
                          END INTERFACE Compute_Solutions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           Contains
!***********************************************************************
!***********************************************************************
!deck read_matrix
!***begin prologue     read_matrix
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            read input or B-spline matrices
!***
!***references
!***routines called
!***end prologue       read_matrix
!
  SUBROUTINE read_matrix
  IMPLICIT NONE
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  LOGICAL                               :: dollar, logkey
  INTEGER                               :: intkey
  INTEGER                               :: i
  REAL*8                                :: fpkey
  INTEGER                               :: len, lenth
  INTEGER                               :: iostat
!
!
  IF ( dollar('$matrix_parameters',card,cpass,inp) ) then
       i0stat=chrkey(card,'initial_state','from_input',' ')
       matrix_file=chrkey(card,'matrix_file_name','from_file',' ')
       matrix_type=chrkey(card,'matrix_type','real',' ')
       triangle=logkey(card,'triangular_storage',.false.,' ')
       packed=logkey(card,'pack_matrices',.false.,' ')
       in_core=logkey(card,'in_core',.false.,' ')
       print_parameter=logkey(card,'print=on',.false.,' ')
       drop_tol = fpkey(card,'drop_tolerance',1.d-08,' ')
       solve_only=logkey(card,'solve_only',.false.,' ')
       rows_to_print = intkey(card,'rows_to_print',20,' ')
       lenbuf = intkey(card,'buffer_length',lenbuf,' ')
       preconditioner=chrkey(card,'preconditioner','none',' ')
       read_only = logkey(card,'read_only',.false.,' ')
       len=lenth(matrix_file)
       WRITE(iout,*) 'Matrix File Name = ',matrix_file(1:len)
  END IF
  IF( matrix_file == 'from_input') THEN
      write(iout,*) 'Matrix File from Input'
      n3d = intkey(card,'matrix_size',2,' ')
  ELSE 
      write(iout,*) 'Matrix File from Disk'
      OPEN(UNIT=50,FILE=matrix_file(1:len),ACCESS='sequential',        &
           FORM='unformatted',IOSTAT=IOSTAT,STATUS='old')
      READ(50) n3d
  END IF
  IF(i0stat /= 'from_disk') THEN
     write(iout,*) 'Opening File to hold data'
     OPEN(UNIT=20,FILE='solution',ACCESS='sequential',                 &
          FORM='unformatted',IOSTAT=IOSTAT,STATUS='unknown')
  END IF
  IF(in_core) THEN
     lenbuf = n3d*(n3d-1)/2
  ELSE
     lenbuf = min(lenbuf,n3d*(n3d-1)/2)
  END IF
  rows_to_print = min(rows_to_print,n3d)      
  ALLOCATE(rowlab(n3d),collab(n3d))
  DO i=1,n3d
     rowlab(i) = 'Row '//itoc(i)
     collab(i) = 'Col '//itoc(i)
  END DO
  WRITE(iout,1) n3d
  WRITE(iout,*) 'Matrix Type = ',matrix_type(1:len)
  IF (matrix_type == 'real') THEN
      IF (triangle) THEN
!
!         Use Triangular Storage.
!
          ALLOCATE(triangle_matrix_d(n3d*(n3d+1)/2), rhs_d(n3d))
          Call Compute_Solutions(triangle_matrix_d,rhs_d) 
      ELSE
!
!         Use Full Storage
!
          ALLOCATE(matrix_d(n3d,n3d), rhs_d(n3d))
          Call Compute_Solutions(matrix_d,rhs_d) 
      END IF
  ELSE
      IF (triangle) THEN
          ALLOCATE(triangle_matrix_z(n3d*(n3d+1)/2), rhs_z(n3d))
          Call Compute_Solutions(triangle_matrix_z,rhs_z) 
      ELSE
          ALLOCATE(matrix_z(n3d,n3d), rhs_z(n3d))
             Call Compute_Solutions(matrix_z,rhs_z) 
      END IF
  END IF
  DEALLOCATE(rowlab,collab)
1 FORMAT(/,10x,'Matrix Size = ',i8)
2 Format(/,20x,'Cholesky Decomposition of Overlap Matrix')
END SUBROUTINE read_matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!***********************************************************************
!deck Compute_Solutions_d
!***begin prologue     Compute_Solutions_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Compute_Solutions_d
!
  SUBROUTINE Compute_Solutions_d (matrix,rhs)
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:)                :: matrix
  REAL*8, DIMENSION(:)                  :: rhs
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  CHARACTER (LEN=8)                     :: kewrd
  INTEGER                               :: intkey
  LOGICAL                               :: dollar, logkey
  REAL*8                                :: fpkey
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: IOSTAT
  INTEGER                               :: count
  REAL*8, DIMENSION(:),   ALLOCATABLE   :: scr_mat
  REAL*8, DIMENSION(:),   ALLOCATABLE   :: scr_rhs
  INTEGER, DIMENSION(:),  ALLOCATABLE   :: ipvt
!
  IF( matrix_file == 'from_input') THEN
      DO i=1,n3d
!
!        Read a row up to the diagonal
!
         call fparr(card,'matrix_row_'//itoc(i),matrix(1:i,i),i,' ')
         matrix(i,1:i) = matrix(1:i,i)
      END DO
      call fparr(card,'right_hand_side',rhs,n3d,' ')      
  ELSE
      DO i=1,n3d
         READ(50) matrix(i,1:i)
         matrix(1:i,i) = matrix(i,1:i)
      END DO
      rhs(:) = 0.d0
      rhs(1) = 1.d0
      rhs(n3d) = 1.d0
  END IF
  IF(print_parameter) THEN          
     title='Input Matrix'
     call matprt(title,matrix,n3d,n3d,rows_to_print,rows_to_print,   &
                 0,0,rowlab,collab,1,eig,ifeig)    
     title='Input Rhs'
     call matprt(title,rhs,n3d,1,rows_to_print,1,0,0,rowlab,collab, &
                 1,eig,ifeig)
  END IF
  IF ( .not.read_only) THEN 
      ALLOCATE( scr_mat(n3d*(n3d+1)/2), scr_rhs(1:n3d), ipvt(n3d) )
!   
!     Put the matrix in packed lower triagular form by columns
!     a(11), a(21),...a(n,1),a(22), a(32),...a(n,2)  etc.
!
      count = 0
      DO i=1,n3d
         DO j=i,n3d
            count = count + 1
            scr_mat(count) = matrix(j,i) 
         END DO
      END DO
      scr_rhs(:)   = rhs(:) 
!
!
      call dspsv('l',n3d,1,scr_mat,ipvt,scr_rhs,n3d,info)
      title='Solution'
      call matprt(title,scr_rhs,n3d,1,n3d,1,0,0,rowlab,collab, &
                  1,eig,ifeig)
      IF(i0stat /= 'from_disk') THEN
         WRITE(20) scr_rhs(:)
         CLOSE(20)  
      END IF
      DEALLOCATE( scr_mat, scr_rhs, ipvt )
  END IF
!***********************************************************************
  END SUBROUTINE Compute_Solutions_d
!***********************************************************************
!***********************************************************************
!deck Compute_Solutions_z
!***begin prologue     Compute_Solutions_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Compute_Solutions_z
!
  SUBROUTINE Compute_Solutions_z(matrix,rhs)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:,:)                :: matrix 
  COMPLEX*16, DIMENSION(:)                  :: rhs 
  CHARACTER(LEN=80)                         :: chrkey
  CHARACTER (LEN=8)                         :: itoc
  CHARACTER (LEN=8)                         :: kewrd
  LOGICAL                                   :: dollar, logkey
  INTEGER                                   :: intkey
  REAL*8                                    :: fpkey
  INTEGER                                   :: i
  INTEGER                                   :: j
  INTEGER                                   :: count
  INTEGER                                   :: iostat
  REAL*8, DIMENSION(:),       ALLOCATABLE   :: rinput
  COMPLEX*16, DIMENSION(:),   ALLOCATABLE   :: scr_mat
  COMPLEX*16, DIMENSION(:),   ALLOCATABLE   :: scr_rhs
  INTEGER, DIMENSION(:),      ALLOCATABLE   :: ipvt
!
  ALLOCATE(rinput(n3d))
  IF( matrix_file == 'from_input') THEN
      DO i=1,n3d
         call fparr(card,'matrix_row_'//itoc(i),rinput(1:i),i,' ')
         matrix(1:i,i) = rinput(1:i) 
         matrix(i,1:i) = matrix(1:i,i)
      END DO
      call fparr(card,'right_hand_side',rinput,n3d,' ')      
      rhs(:) = rinput(:)
  ELSE 
      DO i=1,n3d
         READ(50) rinput(1:i)
         matrix(1:i,i) = rinput(1:i)
         matrix(i,1:i) = matrix(1:i,i)
      END DO
      rhs(:) = (0.d0,0.d0)
      rhs(1) = (1.d0,0.d0)
  END IF
  DEALLOCATE(rinput)
  IF(print_parameter) THEN          
     title='Input Matrix'
     call prntcmn(title,matrix,n3d,rows_to_print,n3d,n3d,iout,'e')
     title='Input Rhs'
     call prntcmn(title,rhs,n3d,1,n3d,1,iout,'e')
  END IF
  IF(.not.read_only) THEN
      ALLOCATE( scr_mat(n3d*(n3d+1)/2), scr_rhs(1:n3d), ipvt(n3d) )
      DO i=1,n3d
         DO j=i,n3d
            count = count + 1
            scr_mat(count) = matrix(j,i) 
         END DO
      END DO
      scr_rhs(:)   = rhs(:) 
      call chpsv('l',n3d,1,scr_mat,ipvt,scr_rhs,n3d,info)
      title='Solution'
      call prntcmn(title,scr_rhs,n3d,1,n3d,1,iout,'e')
      IF(i0stat /= 'from_disk') THEN
         WRITE(20) scr_rhs(:)
         CLOSE(20)  
      END IF
      DEALLOCATE( scr_mat, scr_rhs, ipvt )
  END IF
! ***********************************************************************
  END SUBROUTINE Compute_Solutions_z
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!***********************************************************************
!deck Compute_Solutions_Tri_d
!***begin prologue     Compute_Solutions_Tri_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Compute_Solutions_Tri_d
!
  SUBROUTINE Compute_Solutions_Tri_d (matrix,rhs)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                  :: matrix
  REAL*8, DIMENSION(:)                  :: rhs
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  INTEGER                               :: intkey
  LOGICAL                               :: dollar, logkey
  REAL*8                                :: fpkey
  REAL*8                                :: dum
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: ij
  INTEGER                               :: count
  REAL*8, DIMENSION(:),   ALLOCATABLE   :: scr_mat
  REAL*8, DIMENSION(:),   ALLOCATABLE   :: scr_rhs
  INTEGER, DIMENSION(:),  ALLOCATABLE   :: ipvt
!
  IF( matrix_file == 'from_input') THEN
      count = 0
      DO i=1,n3d
         call fparr(card,'matrix_row_'//itoc(i),matrix(count+1),i,' ')      
         count = count + i
      END DO
  ELSE
      count = 0
      DO i=1,n3d
         READ(50) ( matrix(j), j=count+1, count + i )      
         count = count + i
      END DO
  END IF
  call fparr(card,'right_hand_side',rhs,n3d,' ')      
  IF(print_parameter) THEN          
     title='Input Matrix'
     count = 0
     DO i=1,n3d
        write(iout,1) i
        write(iout,2) ( matrix(j), j=count+1, count + i )
        count = count + i
     END DO
     title='Input Rhs'
     call matprt(title,rhs,n3d,1,1,1,0,0,rowlab,collab,1,eig,ifeig)      
  END IF
  IF(.not.read_only) THEN
      ALLOCATE( scr_mat(n3d*(n3d+1)/2), scr_rhs(n3d), ipvt(n3d) )
      count = 0
      DO i=1,n3d
         DO j=1,i
            count = count + 1
            ij = j*(j-1)/2 + i
            scr_mat(count) = matrix(ij)
         END DO
      END DO
      scr_rhs(:)   = rhs(:) 
      call dspsv('l',n3d,1,scr_mat,ipvt,scr_rhs,n3d,info )
      title='Solution'
      call matprt(title,scr_rhs,n3d,1,rows_to_print,1,0,0,rowlab,collab, &
                  1,eig,ifeig)
      IF(i0stat /= 'from_disk') THEN
         WRITE(20) scr_rhs(:)
         CLOSE(20)  
      END IF
      DEALLOCATE( scr_mat, scr_rhs, ipvt )
  END IF
1 Format(/,5x,'Row = ',i4,/)
2 Format( (15x,5f10.5) )
!***********************************************************************
  END SUBROUTINE Compute_Solutions_Tri_d
!***********************************************************************
!***********************************************************************
!deck Compute_Solutions_Tri_z
!***begin prologue     Compute_Solutions_Tri_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Compute_Solutions_Tri_z
!
  SUBROUTINE Compute_Solutions_Tri_z(matrix,rhs)
                           USE Iterative_Global
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)                 :: matrix 
  COMPLEX*16, DIMENSION(:)                 :: rhs 
  COMPLEX*16                               :: dum 
  CHARACTER(LEN=80)                        :: chrkey
  CHARACTER (LEN=8)                        :: itoc
  LOGICAL                                  :: dollar, logkey
  INTEGER                                  :: intkey
  REAL*8                                   :: fpkey
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: ij
  INTEGER                                  :: count
  REAL*8, DIMENSION(:),    ALLOCATABLE     :: rinput
  COMPLEX*16, DIMENSION(:), ALLOCATABLE    :: scr_mat
  COMPLEX*16, DIMENSION(:),   ALLOCATABLE  :: scr_rhs
  INTEGER, DIMENSION(:),      ALLOCATABLE  :: ipvt
!
  ALLOCATE(rinput(n3d))
  IF( matrix_file == 'from_input') THEN
      count = 0
      DO i=1,n3d
         call fparr(card,'matrix_row_'//itoc(i),rinput,i,' ')      
         DO j=1,i  
            count = count + 1
            matrix(count) = rinput(j)
         END DO
      END DO
  ELSE
      count = 0
      DO i=1,n3d
         READ(50) ( rinput(j), j=1, i )      
         DO j=1,i
            count = count + 1
            matrix(count) = rinput(j)
         END DO
      END DO
  END IF
  call fparr(card,'right_hand_side'//itoc(i),rinput,n3d,' ')      
  rhs(:) = rinput(:)
  IF(print_parameter) THEN          
     title='Input Matrix'
     count = 0
     DO i=1,n3d
        write(iout,1) i
        write(iout,2) (matrix(j), j=count+1, count + i )
        count = count + i
     END DO       
  END IF
  IF ( .not.read_only) THEN 
      ALLOCATE( scr_mat(n3d*(n3d+1)/2), scr_rhs(n3d), ipvt(n3d) )
      count = 0
      DO i=1,n3d
         DO j=1,i
            count = count + 1
            ij = j*(j-1)/2 + i
            scr_mat(count) = matrix(ij)
         END DO
      END DO
      call chpsv('l',n3d,1,scr_mat,ipvt,scr_rhs,n3d,info )
      title='Solution'
      call prntcmn(title,scr_rhs,n3d,1,n3d,1,iout,'e')  
      IF(i0stat /= 'from_disk') THEN      
         WRITE(20) scr_rhs(:)
         CLOSE(20)  
      END IF
      DEALLOCATE( scr_mat, scr_rhs, ipvt )
  END IF
1  Format(/,5x,'Row = ',i4,/)
2  Format( (15x,6f10.5))
3  Format(/,20x,'Cholesky Decomposition of Overlap Matrix')
4  Format(a80)
!***********************************************************************
  END SUBROUTINE Compute_Solutions_Tri_z
!***********************************************************************
  END  MODULE Matrix_Module
!***********************************************************************
!***********************************************************************
