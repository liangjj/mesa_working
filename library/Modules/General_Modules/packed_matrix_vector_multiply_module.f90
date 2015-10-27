!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! packed_matrix_vector_multiply_module
!**begin prologue     packed_matrix_vector_multiply_module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains the routines needed to multiply the
!***                  Hamiltonian on a vector or to mutiply the exponential
!***                  sector propagators on a vector. Explicit interfaces are 
!***                  used to allow a transparent use of generic subroutines 
!***                  which work for both real and complex vectors.  
!***                  This feature permits a single code to be used for both 
!***                  real and imaginary time propagation.
!***description       See the specific routined.
!**references
!**modules needed     See USE statements below
!**end prologue       packed_matrix_vector_multiply_module
!***********************************************************************
!***********************************************************************
                      MODULE packed_matrix_vector_multiply_module
                         USE dvrprop_global
                         USE dvr_shared
                         USE dvr_global
                         USE Iterative_Global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE packed_symmetric_matrix_on_vector
              MODULE PROCEDURE packed_symmetric_matrix_on_vector_d,                      &
                               packed_symmetric_matrix_on_vector_z,                      &
                               packed_symmetric_matrix_on_vector_zz,                     &
                               packed_symmetric_matrix_on_vector_e_d,                    &  
                               packed_symmetric_matrix_on_vector_e_z,                    &    
                               packed_symmetric_matrix_on_vector_e_zz    
                END INTERFACE packed_symmetric_matrix_on_vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE column_packed_symmetric_matrix_on_vector
              MODULE PROCEDURE column_packed_symmetric_matrix_on_vector_d,               &
                               column_packed_symmetric_matrix_on_vector_z,               &
                               column_packed_symmetric_matrix_on_vector_zz
                END INTERFACE column_packed_symmetric_matrix_on_vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE column_packed_matrix_on_vector
              MODULE PROCEDURE column_packed_matrix_on_vector_d,                         &
                               column_packed_matrix_on_vector_z,                         & 
                               column_packed_matrix_on_vector_zz 
                END INTERFACE column_packed_matrix_on_vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE u_l_tran_column_packed_matrix_u_r_on_vector
              MODULE PROCEDURE u_l_tran_column_packed_matrix_u_r_on_vector_d,            &
                               u_l_tran_column_packed_matrix_u_r_on_vector_z
                END INTERFACE u_l_tran_column_packed_matrix_u_r_on_vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE u_l_tran_column_packed_matrix_u_r
              MODULE PROCEDURE u_l_tran_column_packed_matrix_u_r_d,                      &
                               u_l_tran_column_packed_matrix_u_r_z
                END INTERFACE u_l_tran_column_packed_matrix_u_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE column_packed_lower_triangle_on_vector
              MODULE PROCEDURE column_packed_lower_triangle_on_vector_d,                 &
                               column_packed_lower_triangle_on_vector_z,                 &             
                               column_packed_lower_triangle_on_vector_zz             
                END INTERFACE column_packed_lower_triangle_on_vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    INTERFACE column_packed_upper_triangle_on_vector
              MODULE PROCEDURE column_packed_upper_triangle_on_vector_d,                 &
                               column_packed_upper_triangle_on_vector_z,                 &             
                               column_packed_upper_triangle_on_vector_zz               
                END INTERFACE column_packed_upper_triangle_on_vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck packed_symmetric_matrix_on_vector_d
!***begin prologue     packed_symmetric_matrix_on_vector_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           matrix multiply
!***author             schneider, barry (nsf)
!***source             
!***purpose            Multiply the packed hamiltonian matrix on a vector.
!***                   Only the lower triangle of the matrix is stored.
!***                   The packing is assumed to be random with the two non zero
!***                   indices and non zero element stored..
!***references         
!
!***routines called    
!***end prologue       packed_symmetric_matrix_on_vector_d
  Subroutine packed_symmetric_matrix_on_vector_d (hbuf,ibuf,diag,                   &
                                                  vector_in,vector_out,non_zero)
  IMPLICIT NONE
  REAL*8,  DIMENSION(:)        :: hbuf
  INTEGER, DIMENSION(:,:)      :: ibuf
  REAL*8,  DIMENSION(:)        :: vector_in
  REAL*8,  DIMENSION(:)        :: vector_out
  INTEGER                      :: non_zero
  REAL*8, DIMENSION(:)         :: diag
  INTEGER                      :: i
  INTEGER                      :: ii
  INTEGER                      :: jj
  vector_out(:) = 0.d0
  if(non_zero /= 0) then
     DO i=1,non_zero
        ii=ibuf(1,i)
        jj=ibuf(2,i)
        vector_out(ii) = vector_out(ii) + hbuf(i) * vector_in(jj)
!
!       Since only lower triangle is stored we need to use the matrix element
!       twice.
!
        vector_out(jj) = vector_out(jj) + hbuf(i) * vector_in(ii)
     END DO
  END IF
  vector_out(:) = vector_out(:) + diag(:) * vector_in(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE packed_symmetric_matrix_on_vector_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***begin prologue     packed_symmetric_matrix_on_vector_z
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           matrix multiply
!***author             schneider, barry (nsf)
!***source             
!***purpose            Multiply the packed hamiltonian matrix on a vector.
!***                   Only the lower triangle of the matrix is stored.
!***                   See comments on last subroutine.  Here its only necessary
!***                   to do a complex conjugation.
!***references         
!
!***routines called    
!***end prologue       packed_symmetric_matrix_on_vector_z
  Subroutine packed_symmetric_matrix_on_vector_z(hbuf,ibuf,diag,                    &
                                                 vector_in,vector_out,non_zero)
  IMPLICIT NONE
  REAL*8,      DIMENSION(:)      :: hbuf
  INTEGER, DIMENSION(:,:)        :: ibuf
  COMPLEX*16,  DIMENSION(:)      :: vector_in
  COMPLEX*16,  DIMENSION(:)      :: vector_out
  INTEGER                        :: non_zero
  REAL*8,     DIMENSION(:)       :: diag
  INTEGER                        :: i
  INTEGER                        :: ii
  INTEGER                        :: jj
  vector_out(:) = (0.d0,0.d0)
  if(non_zero /= 0) then
     DO i=1,non_zero
        ii=ibuf(1,i)
        jj=ibuf(2,i)
        vector_out(ii) = vector_out(ii) + hbuf(i) * vector_in(jj)
        vector_out(jj) = vector_out(jj) + hbuf(i) * vector_in(ii)
     END DO
  END IF
  vector_out(:) = vector_out(:) + diag(:) * vector_in(:)
END SUBROUTINE packed_symmetric_matrix_on_vector_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***begin prologue     packed_symmetric_matrix_on_vector_zz
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           matrix multiply
!***author             schneider, barry (nsf)
!***source             
!***purpose            Multiply the packed hamiltonian matrix on a vector.
!***                   Only the lower triangle of the matrix is stored.
!***                   See comments on last subroutine.  Here its only necessary
!***                   to do a complex conjugation.
!***references         
!
!***routines called    
!***end prologue       packed_symmetric_matrix_on_vector_zz
  Subroutine packed_symmetric_matrix_on_vector_zz(hbuf,ibuf,diag,                    &
                                                 vector_in,vector_out,non_zero)
  IMPLICIT NONE
  COMPLEX*16,  DIMENSION(:)      :: hbuf
  INTEGER, DIMENSION(:,:)        :: ibuf
  COMPLEX*16,  DIMENSION(:)      :: vector_in
  COMPLEX*16,  DIMENSION(:)      :: vector_out
  INTEGER                        :: non_zero
  COMPLEX*16, DIMENSION(:)       :: diag
  INTEGER                        :: i
  INTEGER                        :: ii
  INTEGER                        :: jj
  vector_out(:) = (0.d0,0.d0)
  if(non_zero /= 0) then
     DO i=1,non_zero
        ii=ibuf(1,i)
        jj=ibuf(2,i)
        vector_out(ii) = vector_out(ii) + hbuf(i) * vector_in(jj)
!
!       Since only lower is stored note conjugation to treat Hermitian matrix.
!
        vector_out(jj) = vector_out(jj) + conjg(hbuf(i)) * vector_in(ii)
     END DO
  END IF
  vector_out(:) = vector_out(:) + diag(:) * vector_in(:)
END SUBROUTINE packed_symmetric_matrix_on_vector_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck packed_symmetric_matrix_on_vector_e_d
!***begin prologue     packed_symmetric_matrix_on_vector_e_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           polynomials
!***author             schneider, barry (nsf)
!***source             
!***purpose            Multiply the packed hamiltonian matrix on a vector.
!***                   Only the lower triangle of the matrix is stored.
!***                   Storage also in random, two two index style.
!***references         
!
!***routines called    
!***end prologue       packed_symmetric_matrix_on_vector_e_d
  Subroutine packed_symmetric_matrix_on_vector_e_d (hbuf,ibuf,diag,                   &
                                                    vector_in,vector_out,             &
                                                    non_zero,n_vec)
  IMPLICIT NONE
  REAL*8,  DIMENSION(:)        :: hbuf
  INTEGER, DIMENSION(:,:)      :: ibuf
  REAL*8,  DIMENSION(:,:)      :: vector_in
  REAL*8,  DIMENSION(:,:)      :: vector_out
  INTEGER                      :: n_vec
  INTEGER                      :: non_zero
  REAL*8, DIMENSION(:)         :: diag
  INTEGER                      :: i
  INTEGER                      :: ii
  INTEGER                      :: jj
  vector_out(:,:) = 0.d0
  if(non_zero /= 0) then
     DO i=1,non_zero
        ii=ibuf(1,i)
        jj=ibuf(2,i)
        vector_out(ii,1:n_vec) = vector_out(ii,1:nvec)                               &
                                      +                                              &
                              hbuf(i) * vector_in(jj,1:n_vec)
        vector_out(jj,1:n_vec) = vector_out(jj,1:n_vec)                              &
                                      +                                              &
                              hbuf(i) * vector_in(ii,1:n_vec)
     END DO
  END IF
  DO i = 1, n_vec
     vector_out(:,i) = vector_out(:,i) + diag(:) * vector_in(:,i)
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE packed_symmetric_matrix_on_vector_e_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine packed_symmetric_matrix_on_vector_e_z(hbuf,ibuf,diag,                    &
                                                   vector_in,vector_out,              &
                                                   non_zero,n_vec)
  IMPLICIT NONE
  REAL*8,      DIMENSION(:)      :: hbuf
  INTEGER, DIMENSION(:,:)        :: ibuf
  COMPLEX*16,  DIMENSION(:,:)    :: vector_in
  COMPLEX*16,  DIMENSION(:,:)    :: vector_out
  INTEGER                        :: n_vec
  INTEGER                        :: non_zero
  REAL*8,     DIMENSION(:)       :: diag
  INTEGER                        :: i
  INTEGER                        :: ii
  INTEGER                        :: jj
  vector_out(:,1:n_vec) = (0.d0,0.d0)
  if(non_zero /= 0) then
     DO i=1,non_zero
        ii=ibuf(1,i)
        jj=ibuf(2,i)
        vector_out(ii,1:n_vec) = vector_out(ii,1:n_vec)                               &
                                          +                                           &
                                  hbuf(i) * vector_in(jj,1:n_vec)
        vector_out(jj,1:n_vec) = vector_out(jj,1:n_vec)                               &
                                         +                                            &
                                 hbuf(i) * vector_in(ii,1:n_vec)
     END DO
  END IF
  DO i = 1, n_vec
     vector_out(:,i) = vector_out(:,i) + diag(:) * vector_in(:,i)
  END DO
END SUBROUTINE packed_symmetric_matrix_on_vector_e_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine packed_symmetric_matrix_on_vector_e_zz(hbuf,ibuf,diag,                   &
                                                   vector_in,vector_out,              &
                                                   non_zero,n_vec)
  IMPLICIT NONE
  COMPLEX*16,  DIMENSION(:)      :: hbuf
  INTEGER, DIMENSION(:,:)        :: ibuf
  COMPLEX*16,  DIMENSION(:,:)    :: vector_in
  COMPLEX*16,  DIMENSION(:,:)    :: vector_out
  INTEGER                        :: n_vec
  INTEGER                        :: non_zero
  COMPLEX*16, DIMENSION(:)       :: diag
  INTEGER                        :: i
  INTEGER                        :: ii
  INTEGER                        :: jj
  vector_out(:,1:n_vec) = (0.d0,0.d0)
  if(non_zero /= 0) then
     DO i=1,non_zero
        ii=ibuf(1,i)
        jj=ibuf(2,i)
        vector_out(ii,1:n_vec) = vector_out(ii,1:n_vec)                               &
                                          +                                           &
                                  hbuf(i) * vector_in(jj,1:n_vec)
        vector_out(jj,1:n_vec) = vector_out(jj,1:n_vec)                               &
                                         +                                            &
                                 conjg(hbuf(i)) * vector_in(ii,1:n_vec)
     END DO
  END IF
  DO i = 1, n_vec
     vector_out(:,i) = vector_out(:,i) + diag(:) * vector_in(:,i)
  END DO
END SUBROUTINE packed_symmetric_matrix_on_vector_e_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck column_packed_matrix_on_vector_d
!***begin prologue     column_packed_matrix_on_vector_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           matrix multiply
!***author             schneider, barry (nsf)
!***source             
!***purpose            Multiply a column packed matrix on a vector.
!***                   All of the non-zero elements in each
!***                   column of the matrix are stored.  How these are
!***                   processed depends on whether you are multiplying
!***                   by the matrix or its transpose.
!
!***routines called    
!***end prologue       column_packed_matrix_on_vector_d
  Subroutine column_packed_matrix_on_vector_d (                             &
                                               vector_in,                   &
                                               vector_out,                  &
                                               packed_columns,              &
                                               non_zero_columns,            &
                                               row_index,                   &
                                               transpose )
  IMPLICIT NONE
  REAL*8, DIMENSION(:)               :: packed_columns
  INTEGER, DIMENSION(:)              :: non_zero_columns
  INTEGER, DIMENSION(:)              :: row_index
  REAL*8,  DIMENSION(:)              :: vector_in
  REAL*8,  DIMENSION(:)              :: vector_out
  REAL*8                             :: mat_el
  INTEGER                            :: i
  INTEGER                            :: j
  INTEGER                            :: ij
  INTEGER                            :: count
  INTEGER                            :: n_vec
  INTEGER                            :: n_col
  LOGICAL                            :: transpose
  n_vec=size(vector_in,1)
  n_col=size(non_zero_columns,1)
  vector_out(1:n_vec) = 0.d0
  count = 0
  IF (transpose == .false. ) THEN
      DO j=1,n_col
         DO i=1,non_zero_columns(j)  
            count = count + 1
            mat_el = packed_columns(count)
            ij = row_index(count)
            vector_out(ij) = vector_out(ij) + mat_el * vector_in(j)
         END DO
      END DO
  ELSE
      DO j=1,n_col
         DO i=1,non_zero_columns(j)  
            count = count + 1
            mat_el = packed_columns(count)
            ij = row_index(count)
            vector_out(j) = vector_out(j) + mat_el * vector_in(ij)
         END DO
      END DO
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE column_packed_matrix_on_vector_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine column_packed_matrix_on_vector_z (                                  &
                                                vector_in,                       &
                                                vector_out,                      &
                                                packed_columns,                  &
                                                non_zero_columns,                &
                                                row_index,                       &
                                                transpose )
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                   :: packed_columns
  INTEGER, DIMENSION(:)                  :: non_zero_columns
  INTEGER, DIMENSION(:)                  :: row_index
  COMPLEX*16,  DIMENSION(:)              :: vector_in
  COMPLEX*16,  DIMENSION(:)              :: vector_out
  REAL*8                                 :: mat_el
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: ij
  INTEGER                                :: count
  INTEGER                                :: n_vec
  INTEGER                                :: n_col
  LOGICAL                                :: transpose
  INTEGER                                :: n
  n_vec=size(vector_in,1)
  n_col=size(non_zero_columns,1)
  vector_out(1:n_vec) = 0.d0
  count = 0
  IF (transpose == .false. ) THEN
      DO j=1,n_col
         DO i=1,non_zero_columns(j)  
            count = count + 1
            mat_el = packed_columns(count)
            ij = row_index(count)
            vector_out(ij) = vector_out(ij) + mat_el * vector_in(j)
         END DO
      END DO
  ELSE
      DO j=1,n_col
         DO i=1,non_zero_columns(j)  
            count = count + 1
            mat_el = packed_columns(count)
            ij = row_index(count)
            vector_out(j) = vector_out(j) + mat_el * vector_in(ij)
         END DO
      END DO
  END IF
END SUBROUTINE column_packed_matrix_on_vector_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine column_packed_matrix_on_vector_zz (                                 &
                                                vector_in,                       &
                                                vector_out,                      &
                                                packed_columns,                  &
                                                non_zero_columns,                &
                                                row_index,                       &
                                                transpose )
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)               :: packed_columns
  INTEGER, DIMENSION(:)                  :: non_zero_columns
  INTEGER, DIMENSION(:)                  :: row_index
  COMPLEX*16,  DIMENSION(:)              :: vector_in
  COMPLEX*16,  DIMENSION(:)              :: vector_out
  COMPLEX*16                             :: mat_el
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: ij
  INTEGER                                :: count
  INTEGER                                :: n_vec
  INTEGER                                :: n_col
  LOGICAL                                :: transpose
  INTEGER                                :: n
  n_vec=size(vector_in,1)
  n_col=size(non_zero_columns,1)
  vector_out(1:n_vec) = 0.d0
  count = 0
  IF (transpose == .false. ) THEN
      DO j=1,n_col
         DO i=1,non_zero_columns(j)  
            count = count + 1
            mat_el = packed_columns(count)
            ij = row_index(count)
            vector_out(ij) = vector_out(ij) + mat_el * vector_in(j)
         END DO
      END DO
  ELSE
      DO j=1,n_col
         DO i=1,non_zero_columns(j)  
            count = count + 1
            mat_el = packed_columns(count)
            ij = row_index(count)
            vector_out(j) = vector_out(j) + mat_el * vector_in(ij)
         END DO
      END DO
  END IF
END SUBROUTINE column_packed_matrix_on_vector_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck column_packed_symmetric_matrix_on_vector_d
!***begin prologue     column_packed_symmetric_matrix_on_vector_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           matrix multiply
!***author             schneider, barry (nsf)
!***source             
!***purpose            Multiply the packed hamiltonian matrix on a vector.
!***                   In the packing only the lower triangle is stored
!***references         
!
!***routines called    
!***end prologue       column_packed_symmetric_matrix_on_vector_d
  Subroutine column_packed_symmetric_matrix_on_vector_d (vector_in,                  &
                                                        vector_out,                  &
                                                        packed_columns,              &
                                                        non_zero_columns,            &
                                                        row_index,                   &
                                                        matrix_diagonal)
  IMPLICIT NONE
  REAL*8, OPTIONAL, DIMENSION(:)     :: matrix_diagonal
  REAL*8, DIMENSION(:)               :: packed_columns
  INTEGER, DIMENSION(:)              :: non_zero_columns
  INTEGER, DIMENSION(:)              :: row_index
  REAL*8,  DIMENSION(:)              :: vector_in
  REAL*8,  DIMENSION(:)              :: vector_out
  REAL*8                             :: mat_el
  INTEGER                            :: i
  INTEGER                            :: j
  INTEGER                            :: ij
  INTEGER                            :: count
  INTEGER                            :: n
  n=size(vector_in,1)
  vector_out(:) = 0.d0
  IF ( PRESENT ( matrix_diagonal ) ) THEN
       vector_out(:) = vector_out(:) + matrix_diagonal(:) * vector_in(:)
  END IF
  count = 0
  DO j=1,n
     DO i=1,non_zero_columns(j)  
        count = count + 1
        mat_el = packed_columns(count)
        ij = row_index(count)
!
!       Since only the triangle of the symmetric matrix is stored
!       we need to use the matrix element twice, for a(i,j) and a(j,i).
!
        vector_out(ij) = vector_out(ij) + mat_el * vector_in(j)
        vector_out(j) = vector_out(j)   + mat_el * vector_in(ij)
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE column_packed_symmetric_matrix_on_vector_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine column_packed_symmetric_matrix_on_vector_z (vector_in,                       &
                                                         vector_out,                      &
                                                         packed_columns,                  &
                                                         non_zero_columns,                &
                                                         row_index,                       &
                                                         matrix_diagonal )
  IMPLICIT NONE
  REAL*8,     OPTIONAL, DIMENSION(:)     :: matrix_diagonal
  REAL*8,     DIMENSION(:)               :: packed_columns
  INTEGER, DIMENSION(:)                  :: non_zero_columns
  INTEGER, DIMENSION(:)                  :: row_index
  COMPLEX*16,  DIMENSION(:)              :: vector_in
  COMPLEX*16,  DIMENSION(:)              :: vector_out
  COMPLEX*16                             :: mat_el
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: ij
  INTEGER                                :: count
  INTEGER                                :: n
  n=size(vector_in,1)
  vector_out(:) = 0.d0
  IF ( PRESENT ( matrix_diagonal ) ) THEN
       vector_out(:) = vector_out(:) + matrix_diagonal(:) * vector_in(:)
  END IF
  count = 0
  DO j=1,n
     DO i=1,non_zero_columns(j)  
        count = count + 1
        mat_el = packed_columns(count)
!
!       Since only the triangle of the symmetric matrix is stored
!       we need to use the matrix element twice, for a(i,j) and a(j,i).
!
        ij = row_index(count)
        vector_out(ij) = vector_out(ij) + mat_el * vector_in(j)
        vector_out(j) = vector_out(j)   + mat_el * vector_in(ij)
     END DO
  END DO
END SUBROUTINE column_packed_symmetric_matrix_on_vector_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine column_packed_symmetric_matrix_on_vector_zz (vector_in,                      &
                                                         vector_out,                      &
                                                         packed_columns,                  &
                                                         non_zero_columns,                &
                                                         row_index,                       &
                                                         matrix_diagonal )
  IMPLICIT NONE
  COMPLEX*16, OPTIONAL, DIMENSION(:)     :: matrix_diagonal
  COMPLEX*16, DIMENSION(:)               :: packed_columns
  INTEGER, DIMENSION(:)                  :: non_zero_columns
  INTEGER, DIMENSION(:)                  :: row_index
  COMPLEX*16,  DIMENSION(:)              :: vector_in
  COMPLEX*16,  DIMENSION(:)              :: vector_out
  COMPLEX*16                             :: mat_el
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: ij
  INTEGER                                :: count
  INTEGER                                :: n
  n=size(vector_in,1)
  vector_out(:) = 0.d0
  IF ( PRESENT ( matrix_diagonal ) ) THEN
       vector_out(:) = vector_out(:) + matrix_diagonal(:) * vector_in(:)
  END IF
  count = 0
  DO j=1,n
     DO i=1,non_zero_columns(j)  
        count = count + 1
        mat_el = packed_columns(count)
!
!       Since only the triangle of the symmetric matrix is stored
!       we need to use the matrix element twice, for a(i,j) and a(j,i).
!
        ij = row_index(count)
        vector_out(ij) = vector_out(ij) + mat_el * vector_in(j)
        vector_out(j) = vector_out(j)   + conjg(mat_el) * vector_in(ij)
     END DO
  END DO
END SUBROUTINE column_packed_symmetric_matrix_on_vector_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck u_l_tran_column_packed_matrix_u_r_on_vector_d
!***begin prologue     u_l_tran_column_packed_matrix_u_r_on_vector_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           polynomials
!***author             schneider, barry (nsf)
!***source             
!***purpose            multiply the packed hamiltonian matrix on a vector.
!***references         
!
!***routines called    
!***end prologue       u_l_tran_column_packed_matrix_u_r_on_vector_d
  Subroutine u_l_tran_column_packed_matrix_u_r_on_vector_d (                          &
                                                             u_l,                     &
                                                             u_r,                     &
                                                             vector_in,               &
                                                             vector_out,              &
                                                             packed_columns,          &
                                                             non_zero_columns,        &
                                                             row_index,               &
                                                             matrix_diagonal)
  IMPLICIT NONE
  REAL*8, OPTIONAL,DIMENSION(:)      :: matrix_diagonal
  REAL*8, DIMENSION(:,:)             :: u_l
  REAL*8, DIMENSION(:,:)             :: u_r
  REAL*8, DIMENSION(:)               :: packed_columns
  INTEGER, DIMENSION(:)              :: non_zero_columns
  INTEGER, DIMENSION(:)              :: row_index
  REAL*8,  DIMENSION(:)              :: vector_in
  REAL*8,  DIMENSION(:)              :: vector_out
  REAL*8,  DIMENSION(:), ALLOCATABLE :: vector_scr
  REAL*8                             :: mat_el
  INTEGER                            :: i
  INTEGER                            :: j
  INTEGER                            :: ij
  INTEGER                            :: count
  INTEGER                            :: l_col
  INTEGER                            :: l_row
  INTEGER                            :: r_row
  INTEGER                            :: r_col
  r_col=size(u_r,2)
  r_row=size(u_r,1)
  l_col=size(u_l,2)
  l_row=size(u_l,1)
  ALLOCATE(vector_scr(1:r_row))
  Call ebc(vector_scr,u_r,vector_in,r_row,r_col,1)
  vector_out(1:l_row) = 0.d0
  count = 0
  DO j=1,r_row
     DO i=1,non_zero_columns(j)  
        count = count + 1
        mat_el = packed_columns(count)
        ij = row_index(count)
        vector_out(ij) = vector_out(ij) + mat_el * vector_scr(j)
     END DO
  END DO
  IF ( PRESENT(matrix_diagonal) ) THEN 
       vector_out(:) = vector_out(:) + matrix_diagonal(:) * vector_scr(:)
  END IF
  vector_scr(1:l_row) = vector_out(1:l_row)
  Call ebtc(vector_out,u_l,vector_scr,l_col,l_row,1)
  DEALLOCATE(vector_scr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE u_l_tran_column_packed_matrix_u_r_on_vector_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck u_l_tran_column_packed_matrix_u_r_on_vector_z
!***begin prologue     u_l_tran_column_packed_matrix_u_r_on_vector_z
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           polynomials
!***author             schneider, barry (nsf)
!***source             
!***purpose            multiply the packed hamiltonian matrix on a vector.
!***references         
!
!***routines called    
!***end prologue       u_l_tran_column_packed_matrix_u_r_on_vector_z
 Subroutine u_l_tran_column_packed_matrix_u_r_on_vector_z (                          &
                                                            u_l,                     &
                                                            u_r,                     &
                                                            vector_in,               &
                                                            vector_out,              &
                                                            packed_columns,          &
                                                            non_zero_columns,        &
                                                            row_index,               &
                                                            matrix_diagonal)
  IMPLICIT NONE
  COMPLEX*16, OPTIONAL,DIMENSION(:)      :: matrix_diagonal
  COMPLEX*16, DIMENSION(:,:)             :: u_l
  COMPLEX*16, DIMENSION(:,:)             :: u_r
  COMPLEX*16, DIMENSION(:)               :: packed_columns
  INTEGER, DIMENSION(:)                  :: non_zero_columns
  INTEGER, DIMENSION(:)                  :: row_index
  COMPLEX*16,  DIMENSION(:)              :: vector_in
  COMPLEX*16,  DIMENSION(:)              :: vector_out
  COMPLEX*16,  DIMENSION(:), ALLOCATABLE :: vector_scr
  COMPLEX*16                             :: mat_el
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: ij
  INTEGER                                :: count
  INTEGER                                :: r_col
  INTEGER                                :: r_row
  INTEGER                                :: l_row
  INTEGER                                :: l_col
  r_col=size(u_r,2)
  r_row=size(u_r,1)
  l_col=size(u_l,2)
  l_row=size(u_l,1)
  ALLOCATE(vector_scr(1:r_row))
  Call cebc(vector_scr,u_r,vector_in,r_row,r_col,1)
  vector_out(1:l_row) = 0.d0
  count = 0
  DO j=1,r_row
     DO i=1,non_zero_columns(j)  
        count = count + 1
        mat_el = packed_columns(count)
        ij = row_index(count)
        vector_out(ij) = vector_out(ij) + mat_el * vector_scr(j)
     END DO
  END DO
  IF ( PRESENT(matrix_diagonal) ) THEN 
       vector_out(:) = vector_out(:) + matrix_diagonal(:) * vector_scr(:)
  END IF
  vector_scr(1:l_row) = vector_out(1:l_row)
  Call cebtc(vector_out,u_l,vector_scr,l_row,l_col,1)
  DEALLOCATE(vector_scr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE u_l_tran_column_packed_matrix_u_r_on_vector_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck u_l_tran_column_packed_matrix_u_r_d
!***begin prologue     u_l_tran_column_packed_matrix_u_r_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           polynomials
!***author             schneider, barry (nsf)
!***source             
!***purpose            multiply the packed hamiltonian matrix on a vector.
!***references         
!
!***routines called    
!***end prologue       u_l_tran_column_packed_matrix_u_r_d
  Subroutine u_l_tran_column_packed_matrix_u_r_d (                                    &
                                                             u_l,                     &
                                                             u_r,                     &
                                                             packed_columns,          &
                                                             non_zero_columns,        &
                                                             row_index,               &
                                                             matrix_out,              &
                                                             matrix_diagonal)
  IMPLICIT NONE
  REAL*8, OPTIONAL,DIMENSION(:)         :: matrix_diagonal
  REAL*8, DIMENSION(:,:)                :: u_l
  REAL*8, DIMENSION(:,:)                :: u_r
  REAL*8, DIMENSION(:)                  :: packed_columns
  INTEGER, DIMENSION(:)                 :: non_zero_columns
  INTEGER, DIMENSION(:)                 :: row_index
  REAL*8,  DIMENSION(:,:)               :: matrix_out
  REAL*8, DIMENSION(:,:), ALLOCATABLE   :: matrix_scratch
  REAL*8                                :: mat_el
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: ij
  INTEGER                               :: count
  INTEGER                               :: l_col
  INTEGER                               :: l_row
  INTEGER                               :: r_row
  INTEGER                               :: r_col
  r_col=size(u_r,2)
  r_row=size(u_r,1)
  l_col=size(u_l,2)
  l_row=size(u_l,1)
  ALLOCATE( matrix_scratch(1:l_row,1:r_col) )
  matrix_scratch(:,:) = 0.d0
  count = 0
  DO j = 1 , r_row
     DO i=1,non_zero_columns(j)  
        count = count + 1
        mat_el = packed_columns(count)
        ij = row_index(count)
        matrix_scratch(ij,1:r_col) = matrix_scratch(ij,1:r_col) + mat_el * u_r(j,1:r_col)
     END DO
  END DO
  IF ( PRESENT(matrix_diagonal) ) THEN 
       DO i = 1, r_col
          matrix_scratch(:,i) = matrix_scratch(:,i) + matrix_diagonal(:) * u_r(:,i)
       END DO
  END IF
  Call ebtc(matrix_out,u_l,matrix_scratch,l_col,l_row,r_col)
  DEALLOCATE( matrix_scratch )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE u_l_tran_column_packed_matrix_u_r_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck u_l_tran_column_packed_matrix_u_r_z
!***begin prologue     u_l_tran_column_packed_matrix_u_r_z
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           polynomials
!***author             schneider, barry (nsf)
!***source             
!***purpose            multiply the packed hamiltonian matrix on a vector.
!***references         
!
!***routines called    
!***end prologue       u_l_tran_column_packed_matrix_u_r_z
 Subroutine u_l_tran_column_packed_matrix_u_r_z (                                    &
                                                            u_l,                     &
                                                            u_r,                     &
                                                            packed_columns,          &
                                                            non_zero_columns,        &
                                                            row_index,               &
                                                            matrix_out,              &
                                                            matrix_diagonal)
  IMPLICIT NONE
  COMPLEX*16, OPTIONAL, DIMENSION(:)          :: matrix_diagonal
  COMPLEX*16, DIMENSION(:,:)                  :: u_l
  COMPLEX*16, DIMENSION(:,:)                  :: u_r
  COMPLEX*16, DIMENSION(:)                    :: packed_columns
  INTEGER, DIMENSION(:)                       :: non_zero_columns
  INTEGER, DIMENSION(:)                       :: row_index
  COMPLEX*16,  DIMENSION(:,:)                 :: matrix_out
  COMPLEX*16,  DIMENSION(:,:), ALLOCATABLE    :: matrix_scratch
  COMPLEX*16                                  :: mat_el
  INTEGER                                     :: i
  INTEGER                                     :: j
  INTEGER                                     :: ij
  INTEGER                                     :: count
  INTEGER                                     :: r_col
  INTEGER                                     :: r_row
  INTEGER                                     :: l_row
  INTEGER                                     :: l_col
  r_col=size(u_r,2)
  r_row=size(u_r,1)
  l_col=size(u_l,2)
  l_row=size(u_l,1)
  ALLOCATE( matrix_scratch(1:l_row,1:r_col) )
  matrix_scratch(:,:) = 0.d0
  count = 0
  DO j = 1 , r_row
     DO i=1,non_zero_columns(j)  
        count = count + 1
        mat_el = packed_columns(count)
        ij = row_index(count)
        matrix_scratch(ij,1:r_col) = matrix_scratch(ij,1:r_col) + mat_el * u_r(j,1:r_col)
     END DO
  END DO
  IF ( PRESENT(matrix_diagonal) ) THEN 
       DO i = 1, r_col
          matrix_scratch(:,i) = matrix_scratch(:,i) + matrix_diagonal(:) * u_r(:,i)
       END DO
  END IF
  Call cebtc(matrix_out,u_l,matrix_scratch,l_col,l_row,r_col)
  DEALLOCATE( matrix_scratch )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE u_l_tran_column_packed_matrix_u_r_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck column_packed_lower_triangle_on_vector_d
!***begin prologue     column_packed_lower_triangle_on_vector_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           polynomials
!***author             schneider, barry (nsf)
!***source             
!***purpose            multiply the column packed lower_triangle on a vector.
!***references         
!
!***routines called    
!***end prologue       column_packed_lower_triangle_on_vector_d
  Subroutine column_packed_lower_triangle_on_vector_d (vector_in,vector_out,      &
                                                       packed_columns,            &
                                                       non_zero_columns,          &
                                                       row_index,                 &
                                                       matrix_diagonal )
  IMPLICIT NONE
  REAL*8, OPTIONAL, DIMENSION(:)     :: matrix_diagonal
  REAL*8, DIMENSION(:)               :: packed_columns
  INTEGER, DIMENSION(:)              :: non_zero_columns
  INTEGER, DIMENSION(:)              :: row_index
  REAL*8,  DIMENSION(:)              :: vector_in
  REAL*8,  DIMENSION(:)              :: vector_out
  REAL*8                             :: mat_el
  INTEGER                            :: n
  INTEGER                            :: i
  INTEGER                            :: j
  INTEGER                            :: ij
  INTEGER                            :: count
  n=size(vector_in,1)
  vector_out(:) = 0.d0
  IF ( PRESENT( matrix_diagonal ) ) THEN
       vector_out(:) = vector_out(:) + matrix_diagonal(:) * vector_in(:)
  END IF
  count = 0
  DO j=1,n
     DO i=1,non_zero_columns(j)  
        count = count + 1
        mat_el = packed_columns(count)
        ij = row_index(count)
        vector_out(ij) = vector_out(ij) + mat_el * vector_in(j)
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE column_packed_lower_triangle_on_vector_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck column_packed_lower_triangle_on_vector_z
!***begin prologue     column_packed_lower_triangle_on_vector_z
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           polynomials
!***author             schneider, barry (nsf)
!***source             
!***purpose            multiply the column packed lower_triangle on a vector.
!***references         
!
!***routines called    
!***end prologue       column_packed_lower_triangle_on_vector_z
  Subroutine column_packed_lower_triangle_on_vector_z (vector_in,vector_out,      &
                                                       packed_columns,            &
                                                       non_zero_columns,          &
                                                       row_index,                 &
                                                       matrix_diagonal )
  IMPLICIT NONE
  REAL*8,     OPTIONAL, DIMENSION(:)     :: matrix_diagonal
  REAL*8,     DIMENSION(:)               :: packed_columns
  INTEGER, DIMENSION(:)                  :: non_zero_columns
  INTEGER, DIMENSION(:)                  :: row_index
  COMPLEX*16,  DIMENSION(:)              :: vector_in
  COMPLEX*16,  DIMENSION(:)              :: vector_out
  INTEGER                                :: n
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: ij
  INTEGER                                :: count
  n=size(vector_in,1)
  vector_out(:) = 0.d0
  IF ( PRESENT( matrix_diagonal ) ) THEN
       vector_out(:) = vector_out(:) + matrix_diagonal(:) * vector_in(:)
  END IF
  count = 0
  DO j=1,n
     DO i=1,non_zero_columns(j)  
        count = count + 1
        ij = row_index(count)
        vector_out(ij) = vector_out(ij) + packed_columns(count) * vector_in(j)
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE column_packed_lower_triangle_on_vector_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck column_packed_lower_triangle_on_vector_zz
!***begin prologue     column_packed_lower_triangle_on_vector_zz
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           polynomials
!***author             schneider, barry (nsf)
!***source             
!***purpose            multiply the column packed lower_triangle on a vector.
!***references         
!
!***routines called    
!***end prologue       column_packed_lower_triangle_on_vector_zz
  Subroutine column_packed_lower_triangle_on_vector_zz (vector_in,vector_out,     &
                                                       packed_columns,            &
                                                       non_zero_columns,          &
                                                       row_index,                 &
                                                       matrix_diagonal )
  IMPLICIT NONE
  COMPLEX*16, OPTIONAL, DIMENSION(:)     :: matrix_diagonal
  COMPLEX*16, DIMENSION(:)               :: packed_columns
  INTEGER, DIMENSION(:)                  :: non_zero_columns
  INTEGER, DIMENSION(:)                  :: row_index
  COMPLEX*16,  DIMENSION(:)              :: vector_in
  COMPLEX*16,  DIMENSION(:)              :: vector_out
  INTEGER                                :: n
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: ij
  INTEGER                                :: count
  n=size(vector_in,1)
  vector_out(:) = 0.d0
  IF ( PRESENT( matrix_diagonal ) ) THEN
       vector_out(:) = vector_out(:) + matrix_diagonal(:) * vector_in(:)
  END IF
  count = 0
  DO j=1,n
     DO i=1,non_zero_columns(j)  
        count = count + 1
        ij = row_index(count)
        vector_out(ij) = vector_out(ij) + packed_columns(count) * vector_in(j)
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE column_packed_lower_triangle_on_vector_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck column_packed_upper_triangle_on_vector_d
!***begin prologue     column_packed_upper_triangle_on_vector_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           polynomials
!***author             schneider, barry (nsf)
!***source             
!***purpose            multiply the column packed upper_triangle on a vector.
!***references         
!
!***routines called    
!***end prologue       column_packed_upper_triangle_on_vector_d
  Subroutine column_packed_upper_triangle_on_vector_d (vector_in,                 &
                                                       vector_out,                &
                                                       packed_columns,            &
                                                       non_zero_columns,          &
                                                       row_index,                 &
                                                       non_zero,                  &
                                                       matrix_diagonal )
  IMPLICIT NONE
  REAL*8, OPTIONAL, DIMENSION(:)     :: matrix_diagonal
  REAL*8, DIMENSION(:)               :: packed_columns
  INTEGER, DIMENSION(:)              :: non_zero_columns
  INTEGER, DIMENSION(:)              :: row_index
  REAL*8,  DIMENSION(:)              :: vector_in
  REAL*8,  DIMENSION(:)              :: vector_out
  REAL*8                             :: mat_el
  INTEGER                            :: non_zero
  INTEGER                            :: n
  INTEGER                            :: i
  INTEGER                            :: j
  INTEGER                            :: ij
  INTEGER                            :: count
  n=size(vector_in,1)
  vector_out(:) = 0.d0
  IF ( PRESENT( matrix_diagonal ) ) THEN
       vector_out(:) = vector_out(:) + matrix_diagonal(:) * vector_in(:)
  END IF
  count = non_zero
  DO j = n, 1, -1
     DO i=1,non_zero_columns(j)  
        mat_el = packed_columns(count)
        ij = row_index(count)
        vector_out(ij) = vector_out(ij) + mat_el * vector_in(j)
        count = count - 1
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE column_packed_upper_triangle_on_vector_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck column_packed_upper_triangle_on_vector_z
!***begin prologue     column_packed_upper_triangle_on_vector_z
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           polynomials
!***author             schneider, barry (nsf)
!***source             
!***purpose            multiply the column packed upper_triangle on a vector.
!***references         
!
!***routines called    
!***end prologue       column_packed_upper_triangle_on_vector_z
  Subroutine column_packed_upper_triangle_on_vector_z (vector_in,                 &
                                                       vector_out,                &
                                                       packed_columns,            &
                                                       non_zero_columns,          &
                                                       row_index,                 &
                                                       non_zero,                  &
                                                       matrix_diagonal )
  IMPLICIT NONE
  REAL*8,     OPTIONAL, DIMENSION(:)     :: matrix_diagonal
  REAL*8,     DIMENSION(:)               :: packed_columns
  INTEGER, DIMENSION(:)                  :: non_zero_columns
  INTEGER, DIMENSION(:)                  :: row_index
  COMPLEX*16,  DIMENSION(:)              :: vector_in
  COMPLEX*16,  DIMENSION(:)              :: vector_out
  INTEGER                                :: non_zero
  INTEGER                                :: n
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: ij
  INTEGER                                :: count
  n=size(vector_in,1)
  vector_out(:) = 0.d0
  IF ( PRESENT( matrix_diagonal ) ) THEN
       vector_out(:) = vector_out(:) + matrix_diagonal(:) * vector_in(:)
  END IF
  count = non_zero
  DO j = n, 1, -1
     DO i=1,non_zero_columns(j)  
        ij = row_index(count)
        vector_out(ij) = vector_out(ij) + packed_columns(count) * vector_in(j)
        count = count - 1
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE column_packed_upper_triangle_on_vector_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck column_packed_upper_triangle_on_vector_zz
!***begin prologue     column_packed_upper_triangle_on_vector_zz
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           polynomials
!***author             schneider, barry (nsf)
!***source             
!***purpose            multiply the column packed upper_triangle on a vector.
!***references         
!
!***routines called    
!***end prologue       column_packed_upper_triangle_on_vector_zz
  Subroutine column_packed_upper_triangle_on_vector_zz (vector_in,                &
                                                       vector_out,                &
                                                       packed_columns,            &
                                                       non_zero_columns,          &
                                                       row_index,                 &
                                                       non_zero,                  &
                                                       matrix_diagonal )
  IMPLICIT NONE
  COMPLEX*16, OPTIONAL, DIMENSION(:)     :: matrix_diagonal
  COMPLEX*16, DIMENSION(:)               :: packed_columns
  INTEGER, DIMENSION(:)                  :: non_zero_columns
  INTEGER, DIMENSION(:)                  :: row_index
  COMPLEX*16,  DIMENSION(:)              :: vector_in
  COMPLEX*16,  DIMENSION(:)              :: vector_out
  INTEGER                                :: non_zero
  INTEGER                                :: n
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: ij
  INTEGER                                :: count
  n=size(vector_in,1)
  vector_out(:) = 0.d0
  IF ( PRESENT( matrix_diagonal ) ) THEN
       vector_out(:) = vector_out(:) + matrix_diagonal(:) * vector_in(:)
  END IF
  count = non_zero
  DO j = n, 1, -1
     DO i=1,non_zero_columns(j)  
        ij = row_index(count)
        vector_out(ij) = vector_out(ij) + packed_columns(count) * vector_in(j)
        count = count - 1
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE column_packed_upper_triangle_on_vector_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END       MODULE packed_matrix_vector_multiply_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
