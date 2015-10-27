!deck h_small_d.f
!***begin prologue     h_small_d
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           small davidson matrix
!***author             schneider, barry (nsf)
!***source
!***
!***references
!***routines called
!***end prologue       h_small_d
  SUBROUTINE h_small_d(h_mat_tri,h_mat_work_tri,iter)
  IMPLICIT NONE
  REAL*8, DIMENSION(0:maximum_number_of_davidson_vectors)      :: h_mat_tri 
  REAL*8, DIMENSION(0:maximum_number_of_davidson_vectors)      :: h_mat_work_tri 
  INTEGER                                                      :: first_vector
  INTEGER                                                      :: last_vector
  INTEGER                                                      :: iter
  INTEGER                                                      :: i
  INTEGER                                                      :: j
  INTEGER                                                      :: ii
  INTEGER                                                      :: index
  INTEGER                                                      :: last_vector
  REAL*8                                                       :: sdot
  CHARACTER (LEN=3)                                            :: itoc
  ii  = first_vector * ( first_vector + 1) / 2
  DO i=first_vector, last_vector
     DO j = 0, i
        index = ii + j
        h_mat_tri(index) = ddot(matrix_size,vec_d(:,i),1,h_vectors_d(:,j),1)
     END DO
     ii = ii + i + 1
  END DO
  h_mat_work_tri(0:index) = h_mat_tri(0:index) 
  IF(prnt) THEN
     title='Small Matrix Iteration = '//itoc(iter)
     write(iout,1) title
     ii = 0
     DO i=0, last_vector
        write(iout,2) i
        write(iout,3) (h_mat_work_tri(j), j=ii, ii + i )
        ii = ii + i + 1
     END DO
  END IF
1 Format(a80)
2 Format(/,5x,'Row = ',i4)
3 Format( (15x,5f10.5) )
END SUBROUTINE h_small_d

