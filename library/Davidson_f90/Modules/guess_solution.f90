!deck @(#)guess_solution
  SUBROUTINE guess_solution(guess_vectors,guess_rhs)
!***begin prologue     guess_solution
!***date written       0803083   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           guess
!***author             
!***source             guess_solution
!***purpose            to form a portion of the matrix and extract starting solutions
!***description
!***references
!***routines called    (none)
!***end prologue       guess_solution
  USE Davidson_Module
  IMPLICIT NONE
  REAL*8, DIMENSION(matrix_size,guess_size)   :: guess_vectors
  REAL*8, DIMENSION(matrix_size)              :: guess_rhs
  INTEGER                                     :: i
  INTEGER                                     :: j
  INTEGER                                     :: ii
  INTEGER                                     :: jj
  INTEGER                                     :: ref_walk
  INTEGER                                     :: i_guess
  INTEGER                                     :: j_guess
  INTEGER                                     :: i_tri
  INTEGER                                     :: count
  REAL*8, DIMENSION(:,:), ALLOCATABLE         :: guess_matrix
  REAL*8, DIMENSION(:), ALLOCATABLE           :: eigenvalues
  REAL*8, DIMENSION(:), ALLOCATABLE           :: temp
  REAL*8, DIMENSION(:), ALLOCATABLE           :: work
  REAL*8                                      :: e_guess
  REAL*8                                      :: tmp
  INTEGER, DIMENSION(:), ALLOCATABLE          :: guess_pointers
  INTEGER, DIMENSION(:), ALLOCATABLE          :: ipivot
!
! This is an attempt to pick out the part of the matrix with the largets diagonals
! It will not always be best, but it cannot be worse.
!
  ALLOCATE(temp(matrix_size), guess_pointers(matrix_size))
  temp(1:matrix_size) = diag_d(1:matrix_size)
  DO i = 1, guess_size
     e_guess = 1.d+30
     DO j = 1, matrix_size
        IF(temp(j) <  e_guess) THEN
           ref_walk = j
           e_guess = temp(j)
        END IF
     END DO
     temp(ref_walk) = 1.0d+30
     guess_pointers(ref_walk) = i
  END DO
!
  IF ( type_calculation = 'linear_system') THEN
       i_tri = guess_size * ( guess_size + 1 ) / 2       
       ALLOCATE( h_mat_tri_d(i_tri) )
       h_mat_tri_d = 0.d0
       DO i = 1, non_zero
          i_guess = guess_pointers(ibuf(1,i))
          IF (i_guess <= 0 ) cycle
          j_guess = guess_pointers(ibuf(2,i))
          IF (j_guess <= 0 ) cycle
          ii = max(i_guess,j_guess)
          jj = max(i_guess,j_guess)
          i_tri = ii * ( ii - 1 ) / 2 + jj      
          h_mat_tri_d(i_tri) = h_mat_tri_d(i_tri)  + h_buf_d(i)
       END DO
       DO  i=1,matrix_size
           i_guess = guess_pointers(i)
           IF (i_guess <= 0) cycle
           i_tri = i_guess * ( i_guess + 1 ) / 2
           h_mat_tri_d(i_tri) = diag_d(i)
       END DO
       IF( prdvd(11) ) THEN
           i_tri = 0
           title='guess Matrix'
           write(iout,1) title
           DO i=1,guess_size
              write(iout,2) i
              write(iout,3) (h_mat_tri_d(j), j=i_tri + 1, i_tri + i )
              i_tri = i_tri + i
           END DO
       END IF
       ALLOCATE(ipivot(guess_size))
       DO  i=1,matrix_size 
           i_guess = guess_pointers(i)
           IF (i_guess <= 0) cycle
           guess_rhs(i_guess) = rhs_d(i)
       END DO
       call dspsv('u',guess_size,1,h_mat_tri_d,ipivot,guess_rhs,guess_size,info)
       temp(1:matrix_size) = 0.d0       
       tmp = 0.0D+00
       DO  i = 1, matrix_size
           IF (guess_pointers(i) > 0) THEN
               temp(i) = guess_rhs(guess_pointers(i))
               IF (ABS(temp(i)) > tmp) THEN
                   tmp = ABS(temp(i))
               END IF
           END IF
       END DO
       guess_rhs(1:matrix_size) = temp(1:matrix_size)
       IF (prdvd(12)) THEN
           title = 'guess right hand side'
           CALL prntrm(title,guess_rhs,guess_size,1,guess_size,1,iout)
       END IF
       DEALLOCATE(temp, guess_pointers, h_mat_tri_d, ipivot)
  ELSE IF( type_calculation = 'eigenvalues') THEN
       ALLOCATE( h_mat_d(guess_size,guess_size)) )
       h_mat_d = 0.d0
       DO i = 1, non_zero
          i_guess = guess_pointers(ibuf(1,i))
          IF (i_guess <= 0 ) cycle
          j_guess = guess_pointers(ibuf(2,i))
          IF (j_guess <= 0 ) cycle
          ii = max(i_guess,j_guess)
          jj = max(i_guess,j_guess)
          h_mat_d(ii,jj) = h_mat_d(ii,jj) + h_buf_d(i)
       END DO
       DO  i=1,matrix_size
           i_guess = guess_pointers(i)
           IF (i_guess <= 0) cycle
           h_mat_d(i_guess,i_guess) = diag_d(i)
       END DO
       DO i = 1, guess_size
          DO j = 1, i
             h_mat_d(j,i) = h_mat_d(i,j)
          END DO
       END DO
       IF (prdvd(11)) THEN
           title = 'guess matrix'
           CALL prntrm(title,h_mat_d,guess_size,guess_size,guess_size,guess_size,iout)
        END IF
       ALLOCATE( work(5*guess_size), eigenvalues(guess_size))
       CALL dsyev('v','l',guess_size,h_mat_d,guess_size,eigenvalues,work,5*guess_size,info)
       DO  i = 1, number_of_guess_vectors
           temp(1:matrix_size) = 0.d0       
           tmp = 0.0D+00
           DO  j = 1, matrix_size
               IF (guess_pointers(j) > 0) THEN
                   temp(j) = h_mat_d(guess_pointers(j),i)
                   IF (ABS(temp(j)) > tmp) THEN
                       tmp = ABS(temp(j))
                   END IF
               END IF
           END DO
           guess_vectors(1:matrix_size,i) = temp(1:matrix_size)
       END DO
       IF (prdvd(12)) THEN
           title = 'guess eigenvalues'
           CALL prntrm(title,eigenvalues,guess_size,1,guess_size,1,iout)
           title = 'guess eigenvectors'
           CALL prntrm(title,h_mat_d,guess_size,guess_size,guess_size,guess_size,iout)
       END IF
       DEALLOCATE(temp, guess_pointers, h_mat_d, work)
  END IF
1    FORMAT (/,t15,'results from guess matrix of size',i5,/,  &
    ' guess   reference   diagonal energy  ' 'guess energy        c(0)')
2    FORMAT (1X,i3,i10,4X,g20.9,g16.9,f8.4)
1 Format(a80)
2 Format(/,5x,'Row = ',i4)
3 Format( (15x,5f10.5) )
END SUBROUTINE guess_solution
