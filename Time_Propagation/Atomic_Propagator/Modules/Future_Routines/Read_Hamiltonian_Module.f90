!***********************************************************************
                           MODULE Read_Hamiltonian_Module
                           USE dvrprop_global
                           USE CC_Prop_Module
                           USE Pack_Hamiltonian_Module
                           USE Preconditioner_Module
                           USE Iterative_Global
                           USE Pack_Global
!
                           IMPLICIT NONE


  CHARACTER (LEN=80)                       :: ham_file, ham_type
  INTEGER                                  :: n_splines
  INTEGER                                  :: n_chan
  INTEGER                                  :: n_bound
  LOGICAL                                  :: ifeig=.false.
  LOGICAL                                  :: diagonalize_only
!***********************************************************************
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          INTERFACE Compute_Eigenstates
                   MODULE PROCEDURE Compute_Eigenstates_d,                               &
                                    Compute_Eigenstates_z
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
  SUBROUTINE read_ham
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
  IF ( dollar('$hamiltonian_parameters',card,cpass,inp) ) then
       ham_file=chrkey(card,'hamiltonian_file_name','from_file',' ')
       ham_type=chrkey(card,'hamiltonian_type','real',' ')
       non_orth=logkey(card,'non_orthogonal_basis',.false.,' ')
       print_parameter=logkey(card,'print=on',.false.,' ')
       diagonalize_only=logkey(card,'diagonalize_only',.false.,' ')
       rows_to_print = intkey(card,'rows_to_print',20,' ')
       eigenvectors_to_print = intkey(card,'eigenvectors_to_print',5,' ')
       preconditioner=chrkey(card,'preconditioner','none',' ')
       packed=logkey(card,'use_packed_matrix',.false.,' ')
       lenbuf=intkey(card,'buffer_length',lenbuf,' ')
       len=lenth(ham_file)
       WRITE(iout,*) 'Hamiltonian File Name = ',ham_file(1:len)
       WRITE(iout,*) 'Non Orthogonal Basis = ',non_orth
  END IF
  IF(preconditioner /= 'none') THEN
     ALLOCATE( number(2), UNIT_NUMBER(2), UNIT_NAME(2))
     UNIT_NUMBER(1) = 1000
     UNIT_NUMBER(2) = 1001
     UNIT_NAME(1) = 'file1000'
     UNIT_NAME(2) = 'file1001'
  END IF
  IF( ham_file == 'from_input') THEN
      n3d = intkey(card,'matrix_size',2,' ')
  ELSE 
      OPEN(UNIT=50,FILE=ham_file(1:len),ACCESS='sequential',                             &
           FORM='unformatted',IOSTAT=IOSTAT,STATUS='old')
      READ(50) n3d
  END IF
  max_row = max(lenbuf/n3d,n3d)
  rows_to_print = min(rows_to_print,n3d)      
  ALLOCATE(rowlab(n3d),collab(n3d))
  DO i=1,n3d
     rowlab(i) = 'Row '//itoc(i)
     collab(i) = 'Col '//itoc(i)
  END DO
  lwork=10*n3d 
  WRITE(iout,1) n3d
  WRITE(iout,*) 'Hamiltonian Type = ',ham_type(1:len)
  IF (ham_type == 'real') THEN
      ALLOCATE(hamiltonian_d(n3d,n3d))
      IF(non_orth) THEN      
         ALLOCATE(overlap_d(n3d,n3d))
      END IF 
      Call Compute_Eigenstates(hamiltonian_d, overlap_d) 
  ELSE
      ALLOCATE(hamiltonian_z(n3d,n3d))
      IF(non_orth) THEN      
         ALLOCATE(overlap_z(n3d,n3d))
      END IF
      Call Compute_Eigenstates(hamiltonian_z, overlap_z) 
  END IF
  DEALLOCATE(rowlab,collab)
  IF(diagonalize_only) THEN
     return
  END IF
1 FORMAT(/,10x,'Matrix Size = ',i8)
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
  INTEGER                               :: intkey
  LOGICAL                               :: dollar, logkey
  REAL*8                                :: fpkey
  INTEGER                               :: i, buflen
!
  IF( ham_file == 'from_input') THEN
      IF(non_orth) THEN      
         ALLOCATE( h_mat_d(1:n3d,1:n3d), eig(1:n3d), work_d(1:lwork), &
                   s_mat_d(1:n3d,1:n3d) )
         DO i=1,n3d
            call fparr(card,'s_'//itoc(i),eig,i,' ')
            s(i,1:i) = eig(1:i)
            s(1:i,i) = eig(1:i)
         END DO
         s_mat_d(1:n3d,1:n3d) = s(1:n3d,1:n3d) 
         IF(print_parameter) THEN          
            title='Input Overlap Matrix'
            call matprt(title,s_mat_d,n3d,n3d,rows_to_print,rows_to_print,0,0,                &
                        rowlab,collab,1,eig,ifeig)
          END IF
      ELSE
         ALLOCATE( h_mat_d(1:n3d,1:n3d), eig(1:n3d), work_d(1:lwork) )
      END IF
      DO i=1,n3d
         call fparr(card,'h_'//itoc(i),eig,i,' ')
         ham(i,1:i) = eig(1:i)
         ham(1:i,i) = eig(1:i)
      END DO
      h_mat_d(1:n3d,1:n3d) = ham(1:n3d,1:n3d) 
      IF(print_parameter) THEN          
         title='Input Hamiltonian Matrix'
         call matprt(title,h_mat_d,n3d,n3d,rows_to_print,rows_to_print,0,0,                   &
                     rowlab,collab,1,eig,ifeig)
      END IF
  ELSE
      IF(non_orth) THEN      
         ALLOCATE( h_mat_d(1:n3d,1:n3d), eig(1:n3d), work_d(1:lwork), &
                   s_mat_d(1:n3d,1:n3d) )
         DO i=1,n3d
            READ(50) s(i,1:i)
            s(1:i,i) = s(i,1:i)
         END DO
!         s(1:n3d,1:n3d) = 0.d0
!         DO i=1,n3d
!            s(i,i) = 1.d0
!         END DO         
         s_mat_d(1:n3d,1:n3d) = s(1:n3d,1:n3d)
         IF(print_parameter) THEN          
            title='Input Overlap Matrix'
            call matprt(title,s_mat_d,n3d,n3d,rows_to_print,rows_to_print,0,0,                &
                        rowlab,collab,1,eig,ifeig)
         END IF
      ELSE
         ALLOCATE( h_mat_d(1:n3d,1:n3d), eig(1:n3d), work_d(1:lwork) )   
      END IF
      DO i=1,n3d
         READ(50) ham(i,1:i)
         ham(1:i,i) = ham(i,1:i)
      END DO
      h_mat_d(1:n3d,1:n3d) = ham(1:n3d,1:n3d)      
      IF(print_parameter) THEN          
         title='Input Hamiltonian Matrix'
         call matprt(title,h_mat_d,n3d,n3d,rows_to_print,rows_to_print,0,0,              &
                     rowlab,collab,1,eig,ifeig)
      END IF
  END IF
  IF(non_orth) THEN
     IF(preconditioner == 'cholesky_decomposition') THEN
        write(iout,1)
        ALLOCATE(lower_d(n3d,n3d)) 
        IF(packed) THEN
           ALLOCATE(non_zero(n3d,2), max_row(2),matrix_diag_d(n3d,2))       
           IF( n3d*n3d <= lenbuf ) THEN
               lenbuf = n3d*n3d
               ALLOCATE (row_buf(lenbuf,2),matrix_buf_d(lenbuf,2))
               in_core=.true.
           ELSE
               ALLOCATE (row_buf(lenbuf,2),matrix_buf_d(lenbuf,2))
               in_core=.false.
           END IF
           unit_pointer = 1
           call h_pack(s_mat_d,.true.)
        END IF
        call Cholesky(s_mat_d,lower_d,n3d)
        title='Cholesky Decomposition'
        call matprt(title,lower_d,n3d,n3d,n3d,n3d,0,0,rowlab,collab,1,eig,ifeig)
        IF(packed) THEN
           unit_pointer = 2
           call h_pack(s_mat_d,.true.)
        END IF
        s_mat_d(1:n3d,1:n3d) = s(1:n3d,1:n3d)
     END IF
     call dsygv(1,'v','l',n3d,h_mat_d,n3d,s_mat_d,n3d,                                   &
                eig,work_d,lwork,info)
  ELSE
     call dsyev('v','l',n3d,h_mat_d,n3d,eig,work_d,lwork,info)
  END IF
  title='Exact Eigenvalues'
  call matprt(title,eig,n3d,1,rows_to_print,1,0,0,rowlab,collab,1,eig,ifeig)
  IF(print_parameter) THEN
     title='Exact Eigenvectors'
     call matprt(title,h_mat_d,n3d,n3d,n3d,eigenvectors_to_print,0,0,                        &
                 rowlab,collab,0,eig,ifeig)
  END IF
  DEALLOCATE( h_mat_d, eig, work_d )
  IF(non_orth) THEN
     DEALLOCATE( s_mat_d )
  END IF
1  FOrmat(/,20x,'Cholesky Decomposition of Overlap Matrix')
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
                           USE Iterative_Global
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:,:)            :: ham 
  COMPLEX*16, DIMENSION(:,:)            :: s 
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  LOGICAL                               :: dollar, logkey
  INTEGER                               :: intkey
  REAL*8                                :: fpkey
  INTEGER                               :: i, buflen
!
  IF( ham_file == 'from_input') THEN
      IF(non_orth) THEN      
         ALLOCATE( h_mat_z(1:n3d,1:n3d), eig(1:n3d), work_z(1:lwork), &
                   rwork(1:lwork), s_mat_z(1:n3d,1:n3d) )
         DO i=1,n3d
            call fparr(card,'s_'//itoc(i),eig,n3d,' ')
            s(i,1:i) = eig(1:i)
            s(1:i,i) = eig(1:i)
         END DO
         s_mat_z(1:n3d,1:n3d) = s(1:n3d,1:n3d) 
         IF(print_parameter) THEN          
            title='Input Overlap Matrix'
            call prntcmn(title,s,n3d,rows_to_print,n3d,n3d,iout,'e')
         END IF
      ELSE
         ALLOCATE( h_mat_z(1:n3d,1:n3d), eig(1:n3d), work_z(1:lwork),   &
                   rwork(1:lwork) )
      END IF
      DO i=1,n3d
         call fparr(card,'h_'//itoc(i),eig,i,' ')
         ham(i,1:i) = eig(1:i)
         ham(1:i,i) = eig(1:i)
      END DO
      h_mat_z(1:n3d,1:n3d) = ham(1:n3d,1:n3d) 
      IF(print_parameter) THEN          
         title='Input Hamiltonian Matrix'
         call prntcmn(title,ham,n3d,rows_to_print,n3d,n3d,n3d,iout,'e')
      END IF
  ELSE 
      IF(non_orth) THEN      
         ALLOCATE( h_mat_z(1:n3d,1:n3d), eig(1:n3d), work_z(1:lwork), &
                   rwork(1:lwork), s_mat_z(1:n3d,1:n3d) )
         DO i=1,n3d
            READ(50) eig(1:i)
            s(i,1:i) = eig(1:i)
            s(1:i,i) = eig(1:i)
         END DO
         s_mat_z(1:n3d,1:n3d) = s(1:n3d,1:n3d)
         IF(print_parameter) THEN          
            title='Input Overlap Matrix'
            call prntcmn(title,s,n3d,rows_to_print,n3d,n3d,iout,'e')
         END IF
      ELSE
         ALLOCATE( h_mat_z(1:n3d,1:n3d), eig(1:n3d), work_z(1:lwork),   &
                    rwork(1:lwork) )   
      END IF
      DO i=1,n3d
         READ(50) eig(1:i)
         ham(i,1:i) = eig(1:i)
         ham(1:i,i) = eig(1:i)
      END DO
      h_mat_z(1:n3d,1:n3d) = ham(1:n3d,1:n3d)      
      IF(print_parameter) THEN          
         title='Input Hamiltonian Matrix'
         call prntcmn(title,ham,n3d,rows_to_print,n3d,n3d,n3d,iout,'e')
      END IF
  END IF
  IF(non_orth) THEN
     IF(preconditioner == 'cholesky_decomposition') THEN
        write(iout,1)
        ALLOCATE(lower_z(n3d,n3d)) 
        IF(packed) THEN
           ALLOCATE(non_zero(n3d,2), max_row(2),matrix_diag_z(n3d,2))        
           IF( n3d*n3d <= lenbuf ) THEN
               ALLOCATE (row_buf(n3d*n3d,2),matrix_buf_z(n3d*n3d,2))
               in_core=.true.
           ELSE
               ALLOCATE (row_buf(lenbuf,2),matrix_buf_z(lenbuf,2))
            in_core=.false.
           END IF
           unit_pointer = 1
           call h_pack(s_mat_z,.true.)
        END IF
        call Cholesky(s_mat_z,lower_z,n3d)
        title='Cholesky Decomposition'
        call prntcmn(title,lower_z,n3d,n3d,n3d,n3d,iout,'e')
        IF(packed) THEN
           unit_pointer = 2
           call h_pack(s_mat_z,.true.)
        END IF
        s_mat_z(1:n3d,1:n3d) = s(1:n3d,1:n3d)
     END IF
     call zhegv(1,'v','l',n3d,h_mat_z,n3d,s_mat_z,n3d,                                &
                eig,work_z,lwork,rwork,info)
  ELSE
     call zheev('v','l',n3d,h_mat_z,n3d,eig,work_z,lwork,rwork,info)
  END IF
  title='Exact Eigenvalues'
  call matprt(title,eig,n3d,1,rows_to_print,1,0,0,rowlab,collab,1,eig,ifeig)
  IF(print_parameter) THEN
     title='Exact Eigenvectors'
     call prntcmn(title,h_mat_z,n3d,n3d,n3d,n3d,iout,'e')
  END IF
  DEALLOCATE( h_mat_z, eig, work_z, rwork )
  IF(non_orth) THEN
     DEALLOCATE( s_mat_z )
  END IF
1  FOrmat(/,20x,'Cholesky Decomposition of Overlap Matrix')
!***********************************************************************
  END SUBROUTINE Compute_Eigenstates_z
!***********************************************************************
  END  MODULE Read_Hamiltonian_Module
!***********************************************************************
!***********************************************************************
