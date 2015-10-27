!***********************************************************************
                           MODULE Direct_Access_IO_Module
                            USE input_output
                            USE dvrprop_global
                            USE Iterative_Global,    ONLY: in_core, lenbuf
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
!***begin prologue     open_direct_access_file
!***date written       000710   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           one-dim
!***author             schneider, barry (nsf)
!***source             
!***purpose             
!***                   
!***references         
!
!***routines called    
!***end prologue       open_direct_access_file
  Subroutine open_direct_access_file(UNIT_NUMBER,UNIT_NAME,record_length)
  IMPLICIT NONE
  INTEGER                      :: UNIT_NUMBER
  CHARACTER(LEN=*)             :: UNIT_NAME
  INTEGER                      :: IOSTAT
  INTEGER                      :: record_length
  write(iout,1)
  write(iout,2) UNIT_NUMBER, UNIT_NAME, record_length
1 FORMAT(/,10x,'Opening a Direct Access File')
2 FORMAT(/,10x,'Unit Number =', i6, 5x, 'Unit Name =',a24, 5x, 'record length in 4 byte words = ',i6)
END SUBROUTINE open_direct_access_file
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
  Subroutine h_pack_z(ham,UNIT_NUMBER,UNIT_NAME,number,nrec,non_zero,row_buf,matrix_buf,prn)
  IMPLICIT NONE
  COMPLEX*16,  DIMENSION(:,:)  :: ham
  COMPLEX*16,  DIMENSION(:)    :: matrix_buf
  INTEGER,  DIMENSION(:)       :: non_zero
  INTEGER,  DIMENSION(:)       :: row_buf
  REAL*8                       :: zero=1.d-14
  INTEGER                      :: i, j
  INTEGER                      :: UNIT_NUMBER
  INTEGER                      :: number
  INTEGER                      :: nrec
  INTEGER                      :: IOSTAT
  LOGICAL                      :: prn
  CHARACTER(LEN=80)            :: title
  CHARACTER(LEN=*)             :: UNIT_NAME
!
!     pack all non-zero, non-diagonal elements
! 
  non_zero(:)=0
  number=0
  IF(in_core) THEN
     DO i=1,n3d
        DO j=1,i
           IF(abs(ham(i,j)) >= zero) THEN 
              non_zero(i) = non_zero(i) + 1
              number = number + 1
              row_buf(number) = j
              matrix_buf = ham(i,j)
           END IF
        END DO
     END DO
  ELSE
     OPEN(UNIT=UNIT_NUMBER,file=UNIT_NAME,ACCESS='sequential',   &
          FORM='unformatted',STATUS='new',IOSTAT=IOSTAT)
     write(UNIT_NUMBER) number, lenbuf
     write(UNIT_NUMBER) non_zero(1:n3d)
     DO i=1,n3d
        DO j=1,i
           IF(abs(ham(i,j)) >= zero) THEN 
              non_zero(i) = non_zero(i) + 1
              number = number + 1
              IF(number > lenbuf) THEN
                 write(UNIT_NUMBER) number
                 write(UNIT_NUMBER) row_buf(1:lenbuf)
                 write(UNIT_NUMBER) matrix_buf(1:lenbuf)
                 number=1
              END IF
              row_buf(number) = j
              matrix_buf(number) = ham(i,j)
           END IF
        END DO
     END DO
     IF(number/= 0) THEN
        write(UNIT_NUMBER) number
        write(UNIT_NUMBER) row_buf(1:lenbuf)
        write(UNIT_NUMBER) matrix_buf(1:lenbuf)
     END IF
     number = 0
     DO i=1,n3d
        number = number + non_zero(i)
     END DO
     rewind(UNIT_NUMBER)        
  END IF
!
! 
  WRITE(iout,1) n3d*(n3d+1)/2, number
1 FORMAT(/,10x,'Possible Non-Zero Elements =', i6,5x, 'Actual Non-Zero Elements =',i6)
END SUBROUTINE h_pack_z
!***********************************************************************
!***********************************************************************
  END  MODULE Pack_Hamiltonian_Module
!***********************************************************************
!***********************************************************************
