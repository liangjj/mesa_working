!***********************************************************************
                           MODULE Pack_Hamiltonian_Module
                            USE input_output
                            USE dvrprop_global
                            USE Iterative_Global
                            USE Pack_Global
                            IMPLICIT NONE
!***********************************************************************
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          INTERFACE h_pack
                   MODULE PROCEDURE h_pack_d,                               &
                                    h_pack_z
                          END INTERFACE h_pack
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           Contains
!***********************************************************************
!***********************************************************************
!***begin prologue     h_pack_d
!***date written       000710   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           one-dim
!***author             schneider, barry (nsf)
!***source             
!***purpose            pack non zero hamiltonian matrix elements and indices 
!***                   
!***references         
!
!***routines called    
!***end prologue       h_pack_d
  Subroutine h_pack_d(upper,prn)
  IMPLICIT NONE
  REAL*8,   DIMENSION(:)       :: upper
  REAL*8                       :: zero=1.d-14
  INTEGER                      :: i
  INTEGER                      :: j
  INTEGER                      :: k
  INTEGER                      :: ipnt
  INTEGER                      :: num
  INTEGER                      :: count
  INTEGER                      :: wptoin
  INTEGER                      :: IOSTAT
  INTEGER                      :: last
  INTEGER                      :: first_row
  INTEGER                      :: last_row
  LOGICAL                      :: prn
  CHARACTER(LEN=80)            :: title
!
!     pack all non-zero elements
! 
  ipnt = unit_pointer
  count = 0  
  non_zero(:,ipnt) = 0
  num = 0
  IF(in_core) THEN
     write(iout,1)
     DO j=1,n3d
        DO i=1,j-1
           num = num + 1
           IF(abs(upper(num)) >= zero) THEN 
              non_zero(i,ipnt) = non_zero(i,ipnt) + 1
              count = count + 1
              row_buf(count,ipnt) = j
              matrix_buf_d(count,ipnt) = upper(num)
           END IF
           num = num + 1
           matrix_diag_d(j,ipnt) = upper(num)
        END DO
     END DO
     number(ipnt)= num
  ELSE
     write(iout,2)
     Call iosys('open '//UNIT_NAME(ipnt)//' as new',                      &
                   0,0,0,UNIT_NAME(ipnt))
     Call iosys('create integer matrix_buffers on '//UNIT_NAME(ipnt),     &
                 -1,0,0,' ')
     DO i=1,n3d
        matrix_diag_d(i,ipnt) = ham(i,i)
     END DO
     trips = n3d/max_row(ipnt)
     last = n3d - trips * max_row(ipnt)
     left_over = max_row(ipnt)
     IF (last /= 0) THEN
         trips = trips + 1
         left_over = last
     END IF         
     first_row = 0
     last_row = 0
     DO i=1,trips
        first_row = last_row + 1
        last_row = last_row + max_row(ipnt)
        IF(i == trips) THEN
           last_row = last_row + left_over
        END IF
        DO j=first_row, last_row
           DO k=1,j-1
              IF(abs(ham(j,k)) >= zero) THEN 
                 num = num + 1
                 count=count+1
                 non_zero(j,ipnt) = non_zero(j,ipnt) + 1
                 row_buf(count,ipnt) = k
                 matrix_buf_d(count,ipnt) = ham(j,k)
              END IF
           END DO
        END DO
        write(iout,3) i, first_row, last_row, count
        Call iosys('write integer matrix_buffers to '                     &
                    //UNIT_NAME(ipnt)//' without rewinding',1,            &
                      count,0,' ')
        Call iosys('write integer matrix_buffers to '                     &
                    //UNIT_NAME(ipnt)//' without rewinding',              &
                      count,row_buf,0,' ')
        Call iosys('write integer matrix_buffers to '                     &
                    //UNIT_NAME(ipnt)//' without rewinding',              &
                      wptoin(count),matrix_buf_d(1,ipnt),0,' ')
        count=0
     END DO
     number(ipnt) = num
     Call iosys('endfile matrix_buffers on '//UNIT_NAME(ipnt),0,0,0,' ')
     Call iosys('write integer total_number_of non_zero_elements to '     &
                 //UNIT_NAME(ipnt),1,number(ipnt),0,' ')
     Call iosys('write integer non_zero_row_indices to '                  &
                 //UNIT_NAME(ipnt),n3d,non_zero(1,ipnt),0,' ')
     Call iosys('write real diagonal_matrix_elements to '                 &
                 //UNIT_NAME(ipnt),n3d,                                   &
                   matrix_diag_d(1,ipnt),0,' ')
     Call iosys('write integer number_of_trips to '                       &
                 //UNIT_NAME(ipnt),1,trips,0,' ')
     Call iosys('write integer left_over to '                             &
                 //UNIT_NAME(ipnt),1,left_over,0,' ')
     Call iosys('rewind all on '//UNIT_NAME(ipnt)//                       &
                ' read-and-write',0,0,0,' ')
!
! 
  END IF
  WRITE(iout,4) n3d*(n3d+1)/2, number
1 FORMAT(/,10x,'Packed Matrix.  All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Packed Matrix.  All Matrix Elements Held on DIsk')
3 FORMAT(/,10x,'Disk Write = ',i5,2x,'First Row = ',i5,2x,'Last Row = ',i5,2x,'Number of Elements = ',i5)
4 FORMAT(/,10x,'Possible Non-Zero Elements =', i6,5x, 'Actual Non-Zero Elements =',i6)
END SUBROUTINE h_pack_d
!***********************************************************************
!***********************************************************************
!***begin prologue     h_pack_z
!***date written       000710   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           one-dim
!***author             schneider, barry (nsf)
!***source             
!***purpose            pack non zero hamiltonian matrix elements and indices 
!***                   
!***references         
!
!***routines called    
!***end prologue       h_pack_z
  Subroutine h_pack_z(ham,prn)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:,:)   :: ham
  REAL*8                       :: zero=1.d-14
  INTEGER                      :: i
  INTEGER                      :: j
  INTEGER                      :: k
  INTEGER                      :: ipnt
  INTEGER                      :: num
  INTEGER                      :: count
  INTEGER                      :: wptoin
  INTEGER                      :: IOSTAT
  INTEGER                      :: last
  INTEGER                      :: first_row
  INTEGER                      :: last_row
  LOGICAL                      :: prn
  CHARACTER(LEN=80)            :: title
!
!     pack all non-zero elements
! 
  ipnt = unit_pointer
  count=0
  non_zero(:,ipnt) = 0
  num = 0
  IF(in_core) THEN
     write(iout,1)
     DO i=1,n3d
        matrix_diag_z(i,ipnt) = ham(i,i)
        DO j=1,i-1
           IF(abs(ham(i,j)) >= zero) THEN 
              non_zero(i,ipnt) = non_zero(i,ipnt) + 1
              num = num + 1
              row_buf(num,ipnt) = j
              matrix_buf_z(num,ipnt) = ham(i,j)
           END IF
        END DO
     END DO
     number(ipnt) = num
  ELSE
     write(iout,2)
     Call iosys('open '//UNIT_NAME(ipnt)//' as new',0,0,0,UNIT_NAME(ipnt))
     Call iosys('create integer matrix_buffers on '//UNIT_NAME(ipnt),-1,0,0,' ')
     DO i=1,n3d
        matrix_diag_z(i,ipnt) = ham(i,i)
     END DO
     trips = n3d/max_row(ipnt)
     last = n3d - trips * max_row(ipnt)
     left_over = max_row(ipnt)
     IF (last /= 0) THEN
         trips = trips + 1
         left_over = last
     END IF         
     first_row = 0
     last_row = 0
     DO i=1,trips
        first_row = last_row + 1
        last_row = last_row + max_row(ipnt)
        IF(i == trips) THEN
           last_row = last_row + left_over
        END IF
        DO j=first_row, last_row
           DO k=1,j-1
              IF(abs(ham(j,k)) >= zero) THEN 
                 num = num + 1
                 count=count+1
                 non_zero(j,ipnt) = non_zero(j,ipnt) + 1
                 row_buf(count,ipnt) = k
                 matrix_buf_z(count,ipnt) = ham(j,k)
              END IF
           END DO
        END DO
        write(iout,3) i, first_row, last_row, count
        Call iosys('write integer matrix_buffers to '                     &
                    //UNIT_NAME(ipnt)//' without rewinding',1,            &
                      count,0,' ')
        Call iosys('write integer matrix_buffers to '                     &
                    //UNIT_NAME(ipnt)//' without rewinding',              &
                      count,row_buf,0,' ')
        Call iosys('write integer matrix_buffers to '                     &
                    //UNIT_NAME(ipnt)//' without rewinding',              &
                      wptoin(2*count),matrix_buf_z(1,ipnt),0,' ')
        count=0
     END DO
     number(ipnt) = num
     Call iosys('endfile matrix_buffers on '//UNIT_NAME(ipnt),0,0,0,' ')
     Call iosys('write integer total_number_of non_zero_elements to '     &
                 //UNIT_NAME(ipnt),1,number(ipnt),0,' ')
     Call iosys('write integer non_zero_row_indices to '                  &
                 //UNIT_NAME(ipnt),n3d,non_zero(1,ipnt),0,' ')
     Call iosys('write real diagonal_matrix_elements to '                 &
                 //UNIT_NAME(ipnt),2*n3d,                                   &
                   matrix_diag_z(1,ipnt),0,' ')
     Call iosys('write integer number_of_trips to '                       &
                 //UNIT_NAME(ipnt),1,trips,0,' ')
     Call iosys('write integer left_over to '                             &
                 //UNIT_NAME(ipnt),1,left_over,0,' ')
     Call iosys('rewind all on '//UNIT_NAME(ipnt)//                       &
                ' read-and-write',0,0,0,' ')
!
! 
! 
  END IF
  WRITE(iout,3) n3d*(n3d+1)/2, number
1 FORMAT(/,10x,'Packed Matrix.  All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Packed Matrix.  All Matrix Elements Held on DIsk')
3 FORMAT(/,10x,'Disk Write = ',i5,2x,'First Row = ',i5,2x,'Last Row = ',i5,2x,'Number of Elements = ',i5)
4 FORMAT(/,10x,'Possible Non-Zero Elements =', i6,5x, 'Actual Non-Zero Elements =',i6)
END SUBROUTINE h_pack_z
!***********************************************************************
!***********************************************************************
  END  MODULE Pack_Hamiltonian_Module
!***********************************************************************
!***********************************************************************
