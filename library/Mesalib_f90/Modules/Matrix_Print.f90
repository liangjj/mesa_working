!***********************************************************************
                           MODULE Matrix_Print
!
                           IMPLICIT NONE
!
!
                         INTERFACE Print_Matrix
                   MODULE PROCEDURE Print_Matrix_d,                                    &
                                    Print_Matrix_z,                                    &
                                    Print_Triangle_Matrix_d,                           &
                                    Print_Triangle_Matrix_z                                     
                          END INTERFACE Print_Matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           Contains
!***********************************************************************
!***********************************************************************
!deck Print_Matrix_d
  subroutine Print_Matrix_d(a,n,m,iout,frmt,title,collab,rowlab)
  implicit none
  REAL*8, DIMENSION(:,:)                   :: a
  INTEGER                                  :: n
  INTEGER                                  :: m
  INTEGER                                  :: iout
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: rowlab
  INTEGER                                  :: ibeg
  INTEGER                                  :: iend
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: k
  INTEGER                                  :: ii
  IF (PRESENT(title)) THEN
      write(iout,1) title
  END IF
  ibeg=0
  IF (PRESENT(frmt)) THEN
      do i=1,m,5
         iend=min(ibeg+5,m)
         ibeg=ibeg+1
         IF (PRESENT(collab)) THEN
             write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         IF ( PRESENT(rowlab) ) THEN
              IF (frmt == 'f') then
                  do j=1,n
                     write (iout,4) rowlab(j),  (a(j,k),k=ibeg,iend)
                 end do
              else
                  do j=1,n
	             write (iout,5) rowlab(j),  (a(j,k),k=ibeg,iend)
                  end do
              end if
         ELSE
              IF (frmt == 'f') then
                  do j=1,n
                     write (iout,6) (a(j,k),k=ibeg,iend)
                  end do
              else
                  do j=1,n
	             write (iout,7)  (a(j,k),k=ibeg,iend)
                  end do
              end if
         END IF
         ibeg=iend
      end do
  ELSE
      do i=1,m,5
         iend=min(ibeg+5,m)
         ibeg=ibeg+1
         IF (PRESENT(collab)) THEN
             write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         IF ( PRESENT(rowlab) ) THEN
              do j=1,n
                 write (iout,5) rowlab(j),  (a(j,k),k=ibeg,iend)
              end do
         ELSE
              do j=1,n
                 write (iout,7) (a(j,k),k=ibeg,iend)
              end do
         END IF
         ibeg=iend
      end do
  END IF
1 format(/,a80,/)
2 format(15x,5(3x,a8,9x))
3 format(1x,'col',5(5x,i6,9x))
4 format(4x,a6,5f20.12)
5 format(4x,a6,5e20.12)
6 format(4x,5f20.12)
7 format(4x,5e20.12)
end subroutine  Print_Matrix_d
!***********************************************************************
!***********************************************************************
!deck Print_Matrix_z
  subroutine Print_Matrix_z(a,n,m,iout,frmt,title,collab,rowlab)
  implicit none
  COMPLEX*16, DIMENSION(:,:)               :: a
  INTEGER                                  :: n
  INTEGER                                  :: m
  INTEGER                                  :: iout
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: rowlab
  INTEGER                                  :: ibeg
  INTEGER                                  :: iend
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: k
  INTEGER                                  :: ii
  IF (PRESENT(title)) THEN
      write(iout,1) title
  END IF
  ibeg=0
  IF (PRESENT(frmt)) THEN
      do i=1,m,2
         iend=min(ibeg+2,m)
         ibeg=ibeg+1
         IF (PRESENT(collab)) THEN
             write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         IF ( PRESENT(rowlab) ) THEN
              IF (frmt == 'f') then
                  do j=1,n
                     write (iout,4) rowlab(j),  (a(j,k),k=ibeg,iend)
                 end do
              else
                  do j=1,n
	             write (iout,5) rowlab(j),  (a(j,k),k=ibeg,iend)
                  end do
              end if
         ELSE
              IF (frmt == 'f') then
                  do j=1,n
                     write (iout,6) (a(j,k),k=ibeg,iend)
                  end do
              else
                  do j=1,n
	             write (iout,7)  (a(j,k),k=ibeg,iend)
                  end do
              end if
         END IF
         ibeg=iend
      end do
  ELSE
      do i=1,m,2
         iend=min(ibeg+2,m)
         ibeg=ibeg+1
         IF (PRESENT(collab)) THEN
             write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         IF ( PRESENT(rowlab) ) THEN
              do j=1,n
                 write (iout,6) rowlab(j),  (a(j,k),k=ibeg,iend)
              end do
         ELSE
              do j=1,n
                 write (iout,7) (a(j,k),k=ibeg,iend)
              end do
         END IF
         ibeg=iend
      end do
  END IF
1 format(/,a80,/)
2 format(15x,2(15x,a8,18x))
3 format(1x,'col(Re,Im)',4(9x,i6,25x))
4 format(4x,a6,4f20.12)
5 format(4x,a6,4e20.12)
6 format(4x,4f20.12)
7 format(4x,4e20.12)
end subroutine  Print_Matrix_z
!***********************************************************************
!***********************************************************************
!deck Print_Triangle_Matrix_d
  subroutine Print_Triangle_Matrix_d(a,n,iout,frmt,title,collab)
  implicit none
  REAL*8, DIMENSION (:)                    :: a
  INTEGER                                  :: n
  INTEGER                                  :: iout
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  INTEGER                                  :: ibeg
  INTEGER                                  :: iend
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: ii
  IF (PRESENT(title)) THEN
      write(iout,1) title
  END IF
  ibeg = 0
  count = 0
  IF (PRESENT(frmt)) THEN
      DO i = 1, n, 5
         iend = min(ibeg+5,n)
         ibeg = ibeg + 1
         IF (PRESENT(collab)) THEN
             write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         IF (frmt == 'f') THEN
             DO j = ibeg, iend
                write (iout,4) a(count + 1:count + j)
                count=count+j
             END DO
         ELSE
             DO j = ibeg, iend
                write (iout,5) a(count + 1:count + j)
                count=count+j
             END DO
         END IF
         ibeg = iend
      END DO
  ELSE
      DO i = 1, n, 5
         iend = min(ibeg+5,n)
         ibeg = ibeg + 1
         IF (PRESENT(collab)) THEN
             write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         DO j = ibeg, iend
            write (iout,5) a(count + 1:count + j)
            count=count+j
         END DO
         ibeg = iend
      END DO
  END IF
1 format(/,a80,/)
2 format(15x,5(3x,a8,9x))
3 format(1x,'col',5(5x,i6,9x))
4 format(10x,5f20.12)
5 format(10x,5e20.12)
end subroutine  Print_Triangle_Matrix_d
!***********************************************************************
!***********************************************************************
!deck Print_Triangle_Matrix_z
  subroutine Print_Triangle_Matrix_z(a,n,iout,frmt,title,collab)
  implicit none
  COMPLEX*16, DIMENSION (:)                :: a
  INTEGER                                  :: n
  INTEGER                                  :: iout
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  INTEGER                                  :: ibeg
  INTEGER                                  :: iend
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: ii
  IF (PRESENT(title)) THEN
      write(iout,1) title
  END IF
  ibeg = 0
  count = 0
  IF (PRESENT(frmt)) THEN
      DO i = 1, n, 5
         iend = min(ibeg+5,n)
         ibeg = ibeg + 1
         IF (PRESENT(collab)) THEN
             write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         IF (frmt == 'f') THEN
             DO j = ibeg, iend
                write (iout,4) a(count + 1:count + j)
                count=count+j
             END DO
         ELSE
             DO j = ibeg, iend
                write (iout,5) a(count + 1:count + j)
                count=count+j
             END DO
         END IF
         ibeg = iend
      END DO
  ELSE
      DO i = 1, n, 5
         iend = min(ibeg+5,n)
         ibeg = ibeg + 1
         IF (PRESENT(collab)) THEN
             write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         DO j = ibeg, iend
            write (iout,5) a(count + 1:count + j)
            count=count+j
         END DO
         ibeg = iend
      END DO
  END IF
1 format(/,a80,/)
2 format(15x,3(3x,a8,9x))
3 format(1x,'col',3(5x,i6,9x))
4 format(10X,5f20.12)
5 format(10x,5e20.12)
end subroutine  Print_Triangle_Matrix_z
!***********************************************************************
!***********************************************************************
END MODULE Matrix_Print
