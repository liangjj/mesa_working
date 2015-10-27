!deck b_small.f
!***begin prologue     b_small
!***date written       980420   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           small davidson matrix
!***author             schneider, barry (nsf)
!***source
!***
!***references
!***routines called
!***end prologue       b_small
  SUBROUTINE b_small(begin,END,action)
INTEGER, INTENT(IN)                      :: begin
INTEGER, INTENT(IN)                      :: END
INTEGER, INTENT(IN OUT)                  :: maxvec
CHARACTER (LEN=*), INTENT(IN)            :: action
IMPLICIT INTEGER (a-z)
REAL*8  sdot
COMMON/io/inp, iout
  IF(action == 'initialize') THEN
     DO i=1,END
        DO j=1,i
           h_mat_d(i,j) = - ddot(matrix_size,vec_d(:,i),1,h_vec_d(:,j),1)
           h_mat_d(j,i) = h_mat_d(i,j)
        END DO
        h_mat_d(i,i) = energy + h_mat_d(i,i)
        DO j=1,number_of_right_hand_sides
           small_rhs_d(i,j) = ddot(matrix_size,vec_d(:,i),1,rhs_d(:,j),1)
        END DO
     END DO
     h_mat_work_d(1:end,1:end) = h_mat_d(1:end,1:end)
     small_rhs_work_d(1:end,1:number_of_right_hand_sides) = small_rhs_d(1:end,1:number_of_right_hand_sides)
  ELSE IF(action == 'fill') THEN
     DO i=1,END
        DO  j=begin,END
            h_mat_d(i,j) = - ddot(matrix_size,vec_d(:,i),1,h_vec_d(:,j),1)
            h_mat_d(j,i) = h_mat_d(i,j)
        END DO
        DO j=1,number_of_right_hand_sides
           small_rhs_d(i,j) = ddot(matrix_size,vec_d(:,i),1,rhs_d(:,j),1)
        END DO
        h_mat_d(i,i) = energy + h_mat_d(i,i)
     END DO
     h_mat_work_d(1:end,begin:end) = h_mat_d(1:end,begin:end)
     small_rhs_work_d(1:end,1:number_of_right_hand_sides) = small_rhs_d(1:end,1:number_of_right_hand_sides)
  END IF
END SUBROUTINE b_small

