!***********************************************************************
                           MODULE Hamiltonian_Module
                           USE dvrprop_global
                           USE CC_Prop_Module
                           USE Preconditioner_Module
                           USE Iterative_Global
                           USE Pack_Global
                           USE Pack_Hamiltonian_Module
!
                           IMPLICIT NONE


  CHARACTER (LEN=80)                       :: ham_file
  INTEGER                                  :: number_of_splines
  INTEGER                                  :: number_of_channels
  INTEGER                                  :: number_of_correlation_terms
  LOGICAL                                  :: ifeig=.false.
  LOGICAL                                  :: diagonalize_only
  LOGICAL                                  :: read_only
!***********************************************************************
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          INTERFACE Compute_Eigenstates
                   MODULE PROCEDURE Compute_Eigenstates_d,                               &
                                    Compute_Eigenstates_Tri_d,                           &
                                    Compute_Eigenstates_z,                               &
                                    Compute_Eigenstates_Tri_z
                          END INTERFACE Compute_Eigenstates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           Contains
!***********************************************************************
!***********************************************************************
!deck read_ham
!***begin prologue     read_ham
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            read input or B-spline Hamiltonian
!***
!***references
!***routines called
!***end prologue       read_ham
!
  SUBROUTINE read_ham(ham_type,type_storage)
  IMPLICIT NONE
  CHARACTER (LEN=*)                     :: ham_type
  CHARACTER (LEN=*)                     :: type_storage
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
  IF ( dollar('$hamiltonian_parameters',card,cpass,inp) ) then
       i0stat=chrkey(card,'initial_state','from_input',' ')
       ham_file=chrkey(card,'hamiltonian_file_name','from_file',' ')
       diagonalize_only=logkey(card,'diagonalize_only',.false.,' ')
       rows_to_print = intkey(card,'rows_to_print',20,' ')
       eigenvectors_to_print = intkey(card,'eigenvectors_to_print',5,' ')
       read_only = logkey(card,'read_only',.false.,' ')
       len=lenth(ham_file)
       WRITE(iout,*) 'Hamiltonian File Name = ',ham_file(1:len)
       WRITE(iout,*) 'Non Orthogonal Basis = ',non_orth
  END IF
  IF( ham_file == 'from_input') THEN
      write(iout,*) 'Hamiltonian File from Input'
      n3d = intkey(card,'matrix_size',2,' ')
  ELSE 
      write(iout,*) 'Hamiltonian File from Disk'
      OPEN(UNIT=50,FILE=ham_file(1:len),ACCESS='sequential',                             &
           FORM='unformatted',IOSTAT=IOSTAT,STATUS='old')
      READ(50) n3d
  END IF
  IF (packed) THEN
      WRITE(iout,1)
  END IF
  IF(i0stat /= 'from_disk') THEN
     write(iout,*) 'Opening File eigenstates to hold data'
     OPEN(UNIT=20,FILE='eigenstates',ACCESS='sequential', FORM='unformatted',            &
          IOSTAT=IOSTAT,STATUS='unknown')
  END IF
  IF(in_core) THEN
     lenbuf = n3d*(n3d-1)/2
  ELSE
     lenbuf = min(lenbuf,n3d*(n3d-1)/2)
  END IF
  rows_to_print = min(rows_to_print,n3d)      
  eigenvectors_to_print = min(eigenvectors_to_print,n3d)
  ALLOCATE(rowlab(n3d),collab(n3d))
  DO i=1,n3d
     rowlab(i) = 'Row '//itoc(i)
     collab(i) = 'Col '//itoc(i)
  END DO
  lwork=10*n3d 
  WRITE(iout,2) n3d
  WRITE(iout,*) 'Hamiltonian Type = ',ham_type(1:len)
  IF (ham_type == 'real') THEN
!
!     
      IF (type_storage == 'triangle') THEN
          ALLOCATE(triangle_hamiltonian_d(n3d*(n3d+1)/2))
          IF(non_orth) THEN      
             ALLOCATE(triangle_overlap_d(n3d*(n3d+1)/2))
          END IF
          Call Compute_Eigenstates(triangle_hamiltonian_d, triangle_overlap_d) 
      ELSE IF(type_storage == 'full') THEN
          ALLOCATE(hamiltonian_d(n3d,n3d))
          IF(non_orth) THEN      
             ALLOCATE(overlap_d(n3d,n3d))
          END IF
          Call Compute_Eigenstates(hamiltonian_d, overlap_d) 
      END IF
  ELSE IF (ham_type == 'complex') THEN
      IF (type_storage == 'triangle') THEN
          ALLOCATE(triangle_hamiltonian_z(n3d*(n3d+1)/2))
          IF(non_orth) THEN      
             ALLOCATE(triangle_overlap_z(n3d*(n3d+1)/2))
          END IF
          Call Compute_Eigenstates(triangle_hamiltonian_z, triangle_overlap_z) 
      ELSE IF(type_storage == 'full') THEN
          ALLOCATE(hamiltonian_z(n3d,n3d))
          IF(non_orth) THEN      
             ALLOCATE(overlap_z(n3d,n3d))
          END IF
          Call Compute_Eigenstates(hamiltonian_z, overlap_z) 
      END IF
  END IF
  DEALLOCATE(rowlab,collab)
  IF(diagonalize_only) THEN
     return
  END IF
1 Format(/,20x,'Using Packed Forms for all Matrices')
2 FORMAT(/,10x,'Matrix Size = ',i8)
3 Format(/,20x,'Cholesky Decomposition of Overlap Matrix')
END SUBROUTINE read_ham
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!***********************************************************************
!deck Compute_Eigenstates_d
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
  REAL*8, DIMENSION(:,:)                :: ham
  REAL*8, DIMENSION(:,:)                :: s
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  CHARACTER (LEN=8)                     :: kewrd
  INTEGER                               :: intkey
  LOGICAL                               :: dollar, logkey
  REAL*8                                :: fpkey
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: ij
  INTEGER                               :: ic
  INTEGER                               :: jc
  INTEGER                               :: is
  INTEGER                               :: js
  INTEGER                               :: first
  INTEGER                               :: last
  INTEGER                               :: count
  INTEGER                               :: iostat
!
  IF( ham_file == 'from_input') THEN
      IF(non_orth) THEN
         DO i=1,n3d
            call fparr(card,'s_'//itoc(i),s(1:i,i),i,' ')
            s(i,1:i) = s(1:i,i)
         END DO
         IF(print_parameter) THEN          
            title='Input Overlap Matrix'
            call matprt(title,s,n3d,n3d,rows_to_print,rows_to_print,0,0,           &
                        rowlab,collab,1,eig,ifeig)
         END IF
      END IF
      DO i=1,n3d
         call fparr(card,'h_'//itoc(i),ham(1:i,i),i,' ')
         ham(i,1:i) = ham(1:i,i)
      END DO
      IF(print_parameter) THEN          
         title='Input Hamiltonian Matrix'
         call matprt(title,ham,n3d,n3d,rows_to_print,rows_to_print,0,0,            &
                     rowlab,collab,1,eig,ifeig)
      END IF
  ELSE
      IF (non_orth) THEN
          DO i=1,n3d
             READ(50) s(i,1:i)
             s(1:i,i) = s(i,1:i)
          END DO
          IF(print_parameter) THEN          
             title='Input Overlap Matrix'
             call matprt(title,s,n3d,n3d,rows_to_print,rows_to_print,0,0,          &
                         rowlab,collab,1,eig,ifeig)
          END IF
      END IF
      DO i=1,n3d
         READ(50) ham(i,1:i)
         ham(1:i,i) = ham(i,1:i)
      END DO
      IF(print_parameter) THEN          
         title='Input Hamiltonian Matrix'
         call matprt(title,ham,n3d,n3d,rows_to_print,rows_to_print,0,0,            &
                     rowlab,collab,1,eig,ifeig)    
      END IF
  END IF
  IF ( .not.read_only) THEN 
      ALLOCATE( h_mat_d(1:n3d,1:n3d), eig(1:n3d), work_d(1:lwork) )
      h_mat_d(1:n3d,1:n3d) = ham(1:n3d,1:n3d) 
      IF(non_orth) THEN      
         ALLOCATE( s_mat_d(1:n3d,1:n3d) )
         s_mat_d(1:n3d,1:n3d) = s(1:n3d,1:n3d)
      END IF
      IF(non_orth) THEN
         call dsygv(1,'v','u',n3d,h_mat_d,n3d,s_mat_d,n3d,eig,work_d,lwork,info)
      ELSE
         call dsyev('v','u',n3d,h_mat_d,n3d,eig,work_d,lwork,info)
      END IF
      title='Exact Eigenvalues'
      call matprt(title,eig,n3d,1,rows_to_print,1,0,0,rowlab,collab,1,eig,ifeig)
      IF(print_parameter) THEN
         title='Exact Eigenvectors'
         call matprt(title,h_mat_d,n3d,n3d,n3d,eigenvectors_to_print,0,0,         &
                     rowlab,collab,0,eig,ifeig) 
      END IF
      IF (i0stat /= 'from_disk') THEN
          kewrd='real'
          write(iout,*) 'Writing '//kewrd//' data to eigenstates'
          write(20) kewrd
          write(20) ( eig(i), i=1,n3d)
          DO i=1,n3d
             write(20) (h_mat_d(j,i), j=1,n3d)
          END DO
          CLOSE(20)
      END IF
      DEALLOCATE( h_mat_d, eig, work_d )
      IF(non_orth) THEN      
         DEALLOCATE( s_mat_d )
      END IF
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
  COMPLEX*16, DIMENSION(:,:)            :: ham 
  COMPLEX*16, DIMENSION(:,:)            :: s 
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  CHARACTER (LEN=8)                     :: kewrd
  LOGICAL                               :: dollar, logkey
  INTEGER                               :: intkey
  REAL*8                                :: fpkey
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: ij
  INTEGER                               :: ic
  INTEGER                               :: jc
  INTEGER                               :: is
  INTEGER                               :: js
  INTEGER                               :: first
  INTEGER                               :: last
  INTEGER                               :: count
  INTEGER                               :: iostat
  REAL*8, DIMENSION(:), ALLOCATABLE     :: rinput
!
  ALLOCATE(rinput(n3d))
  IF( ham_file == 'from_input') THEN
      IF(non_orth) THEN      
         DO i=1,n3d
            call fparr(card,'s_'//itoc(i),rinput(1:i),i,' ')
            s(1:i,i) = rinput(1:i) 
            s(i,1:i) = s(1:i,i)
         END DO
         IF(print_parameter) THEN          
            title='Input Overlap Matrix'
            call prntcmn(title,s,n3d,rows_to_print,n3d,n3d,iout,'e')
         END IF
      END IF
      DO i=1,n3d
         call fparr(card,'h_'//itoc(i),rinput(1:i),i,' ')
         ham(1:i,i) = rinput(1:i) 
         ham(i,1:i) = ham(1:i,i)
      END DO
      IF(print_parameter) THEN          
         title='Input Hamiltonian Matrix'
         call prntcmn(title,ham,n3d,rows_to_print,n3d,n3d,iout,'e')
      END IF
  ELSE 
      IF(non_orth) THEN      
         DO i=1,n3d
            READ(50) rinput(1:i)
            s(1:i,i) = rinput(1:i) 
            s(i,1:i) = s(1:i,i)
         END DO
         IF(print_parameter) THEN          
            title='Input Overlap Matrix'
            call prntcmn(title,s,n3d,rows_to_print,n3d,n3d,iout,'e')
         END IF
      END IF
      DO i=1,n3d
         READ(50) rinput(1:i)
         ham(1:i,i) = rinput(1:i)
         ham(i,1:i) = ham(1:i,i)
      END DO
      IF(print_parameter) THEN          
         title='Input Hamiltonian Matrix'
         call prntcmn(title,ham,n3d,rows_to_print,n3d,n3d,iout,'e')
      END IF
  END IF
  DEALLOCATE(rinput)
  IF(.not.read_only) THEN
      ALLOCATE( h_mat_z(1:n3d,1:n3d), eig(1:n3d), work_z(1:lwork),                     &
                rwork(1:lwork) )
      h_mat_z(1:n3d,1:n3d) = ham(1:n3d,1:n3d)      
      IF(non_orth) THEN      
         ALLOCATE( s_mat_z(1:n3d,1:n3d) )
         s_mat_z(1:n3d,1:n3d) = s(1:n3d,1:n3d) 
      END IF
     IF(non_orth) THEN
        call zhegv(1,'v','u',n3d,h_mat_z,n3d,s_mat_z,n3d,                              &
                   eig,work_z,lwork,rwork,info)
     ELSE
        call zheev('v','u',n3d,h_mat_z,n3d,eig,work_z,lwork,rwork,info)
     END IF
     title='Exact Eigenvalues'
     call matprt(title,eig,n3d,1,rows_to_print,1,0,0,rowlab,collab,1,eig,ifeig)
     IF(print_parameter) THEN
        title='Exact Eigenvectors'
        call prntcmn(title,h_mat_z,n3d,n3d,n3d,n3d,iout,'e')
     END IF
     IF (i0stat /= 'from_disk') THEN
         kewrd='complex'
         write(iout,*) 'Writing '//kewrd//' data to eigenstates'
         write(20) kewrd
         write(20) ( eig(i), i=1,n3d)
         DO i=1,n3d
            write(20) (h_mat_z(j,i), j=1,n3d)
         END DO
         CLOSE(20)
     END IF
     DEALLOCATE( h_mat_z, eig, work_z, rwork )
     IF(non_orth) THEN      
        DEALLOCATE( s_mat_z)
     END IF
  END IF
! ***********************************************************************
  END SUBROUTINE Compute_Eigenstates_z
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!***********************************************************************
!deck Compute_Eigenstates_Tri_d
!***begin prologue     Compute_Eigenstates_Tri_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Compute_Eigenstates_Tri_d
!
  SUBROUTINE Compute_Eigenstates_Tri_d (ham,s)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                  :: ham
  REAL*8, DIMENSION(:)                  :: s
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  INTEGER                               :: intkey
  LOGICAL                               :: dollar, logkey
  REAL*8                                :: fpkey
  REAL*8                                :: dum
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: count
!
  IF( ham_file == 'from_input') THEN
      IF(non_orth) THEN      
         count = 0
         DO i=1,n3d
            call fparr(card,'s_'//itoc(i),s(count+1),i,' ')
            count = count + i   
         END DO
         IF(print_parameter) THEN          
            count = 0
            title='Input Overlap Matrix'
             DO i=1,n3d
                write(iout,1) i
                write(iout,2) (s(j), j=count+1, count + i )
                 count = count + i
             END DO       
          END IF
      END IF
      count = 0
      DO i=1,n3d
         call fparr(card,'h_'//itoc(i),ham(count+1),i,' ')      
         count = count + i
      END DO
  ELSE
      IF(non_orth) THEN      
         count = 0
         DO i=1,n3d
            READ(50) ( s(j), j=count+1, count + i )
            count = count + i
         END DO
         IF(print_parameter) THEN          
            title='Input Overlap Matrix'
            count = 0
            DO i=1,n3d
             write(iout,1) i
             write(iout,2) ( s(j), j=count+1, count + i )
             count = count + i
            END DO
         END IF
      END IF
      count = 0
      DO i=1,n3d
         READ(50) ( ham(j), j=count+1, count + i )      
         count = count + i
      END DO
  END IF
  IF(print_parameter) THEN          
     title='Input Hamiltonian Matrix'
     count = 0
     DO i=1,n3d
        write(iout,1) i
        write(iout,2) (ham(j), j=count+1, count + i )
        count = count + i
     END DO       
  END IF
  IF(.not.read_only) THEN
     ALLOCATE( h_mat_tri_d(n3d*(n3d+1)/2), eig(1:n3d), work_d(1:lwork) )
     h_mat_tri_d(:) = ham(:) 
     IF(non_orth) THEN      
        ALLOCATE( s_mat_tri_d(n3d*(n3d+1)/2) )
        s_mat_tri_d(:) = s(:)
     END IF
     IF(non_orth) THEN
        call dspgv(1,'n','u',n3d,h_mat_tri_d,s_mat_tri_d,eig,dum,n3d,work_d,info)
     ELSE
        call dspev('n','u',n3d,h_mat_tri_d,eig,dum,n3d,work_d,info)
     END IF
     title='Exact Eigenvalues'
     call matprt(title,eig,n3d,1,rows_to_print,1,0,0,rowlab,collab,1,eig,ifeig)  
     DEALLOCATE( h_mat_tri_d, eig, work_d)
     IF(non_orth) THEN      
        DEALLOCATE( s_mat_tri_d )
     END IF
  END IF
1 Format(/,5x,'Row = ',i4,/)
2 Format( (15x,5f10.5) )
!***********************************************************************
  END SUBROUTINE Compute_Eigenstates_Tri_d
!***********************************************************************
!***********************************************************************
!deck Compute_Eigenstates_Tri_z
!***begin prologue     Compute_Eigenstates_Tri_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Compute_Eigenstates_Tri_z
!
  SUBROUTINE Compute_Eigenstates_Tri_z(ham,s)
                           USE Iterative_Global
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)              :: ham 
  COMPLEX*16, DIMENSION(:)              :: s 
  COMPLEX*16                            :: dum 
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  LOGICAL                               :: dollar, logkey
  INTEGER                               :: intkey
  REAL*8                                :: fpkey
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: count
  REAL*8, DIMENSION(:), ALLOCATABLE     :: rinput
!
  ALLOCATE(rinput(n3d))
  IF( ham_file == 'from_input') THEN
      IF(non_orth) THEN      
         count = 0
         DO i=1,n3d
            call fparr(card,'s_'//itoc(i),rinput,i,' ')
            DO j=1,i
               count = count + 1   
               s(count) = rinput(j)
            END DO
         END DO
         IF(print_parameter) THEN          
            title='Input Overlap Matrix'
            count = 0
             DO i=1,n3d
                write(iout,1) i
                write(iout,2) (s(j), j=count+1, count + i )
                count = count + i
             END DO       
         END IF
      END IF
      count = 0
      DO i=1,n3d
         call fparr(card,'h_'//itoc(i),rinput,i,' ')      
         DO j=1,i
            count = count + 1
            ham(count) = rinput(j) 
         END DO
      END DO
  ELSE 
      IF(non_orth) THEN      
         count = 0 
         DO i=1,n3d
            READ(50) rinput(1:i)
            DO j=1,i
              count = count + 1
              s(count) = rinput(j) 
            END DO
         END DO
      END IF
      IF(print_parameter) THEN          
         title='Input Overlap Matrix'
         count = 0
          DO i=1,n3d
             write(iout,1) i
             write(iout,2) (s(j), j=count+1, count + i )
             count = count + i
          END DO       
      END IF
      count = 0
      DO i=1,n3d
         READ(50) rinput(1:i)
         DO j=1,i
            count = count + 1
            ham(count) = rinput(j)
         END DO
      END DO
  END IF
  DEALLOCATE(rinput)
  IF(print_parameter) THEN          
     title='Input Hamiltonian Matrix'
     count = 0
     DO i=1,n3d
        write(iout,1) i
        write(iout,2) (ham(j), j=count+1, count + i )
        count = count + i
      END DO       
  END IF
  IF ( .not.read_only) THEN 
     ALLOCATE( h_mat_tri_z(n3d*(n3d+1)/2), eig(1:n3d), work_z(1:lwork), rwork(1:lwork) )
     h_mat_tri_z(:) = ham(:)      
     IF(non_orth) THEN      
        ALLOCATE( s_mat_tri_z(n3d*(n3d+1)/2) )
        s_mat_tri_z(:) = s(:) 
     END IF
     IF(non_orth) THEN
        call zhpgv(1,'n','u',n3d,h_mat_tri_z,s_mat_tri_z,eig,dum,n3d,work_z,rwork,info)
     ELSE
        call zhpev('n','u',n3d,h_mat_tri_z,eig,dum,n3d,work_z,rwork,info)
     END IF
     title='Exact Eigenvalues'
     call matprt(title,eig,n3d,1,rows_to_print,1,0,0,rowlab,collab,1,eig,ifeig)
     DEALLOCATE( h_mat_tri_z, eig, work_z, rwork )
     IF(non_orth) THEN      
        DEALLOCATE( s_mat_tri_z )
     END IF
  END IF
1  Format(/,5x,'Row = ',i4,/)
2  Format( (15x,6f10.5))
!***********************************************************************
  END SUBROUTINE Compute_Eigenstates_Tri_z
!***********************************************************************
!***********************************************************************
  END  MODULE Hamiltonian_Module
!***********************************************************************
!***********************************************************************
