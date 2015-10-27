!***********************************************************************
                           MODULE Matrix_Print
!
                        USE input_output
                        USE accuracy
                        IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
                     TYPE Real_Matrix
                          CHARACTER(LEN=16)     :: type
                     END TYPE Real_Matrix

                     TYPE Complex_Matrix
                          CHARACTER(LEN=16)     :: type
                     END TYPE Complex_Matrix

                     TYPE Real_Vector
                          CHARACTER(LEN=16)     :: type
                     END TYPE Real_Vector

                     TYPE Complex_Vector
                         CHARACTER(LEN=16)     :: type
                     END TYPE Complex_Vector

                     TYPE Real_Triangle
                         CHARACTER(LEN=16)     :: type
                     END TYPE Real_Triangle

                     TYPE Complex_Triangle
                          CHARACTER(LEN=16)    :: type
                     END TYPE Complex_Triangle
!
!

                         INTERFACE Print_Matrix
                   MODULE PROCEDURE Print_Matrix_d,                                    &
                                    Print_Matrix_z,                                    &
                                    Print_Triangle_Matrix_d,                           &
                                    Print_Triangle_Matrix_z,                           &          
                                    Print_Vector_d,                                    &
                                    Print_Vector_z                                     
                          END INTERFACE Print_Matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           Contains
!***********************************************************************
!***********************************************************************
!deck Print_Matrix_d
  subroutine Print_Matrix_d(type_real_matrix,a,n,m,frmt,title,collab,rowlab)
  implicit none
  TYPE(Real_Matrix)                        :: type_real_matrix
  REAL(idp), DIMENSION(:,:)                :: a
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: rowlab
  INTEGER                                  :: n
  INTEGER                                  :: m
  INTEGER                                  :: lenth
  INTEGER                                  :: ibeg
  INTEGER                                  :: iend
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: k
  INTEGER                                  :: ii
  INTEGER                                  :: len
!
  IF (PRESENT(title)) THEN
      len=lenth(title)
      write(iout,*) title(1:len)
  END IF
  ibeg=0
  IF (PRESENT(frmt)) THEN
      do i=1,m,5
         iend=min(ibeg+5,m)
         ibeg=ibeg+1
         IF (PRESENT(collab)) THEN
             write(iout,1) collab(ibeg:iend)
         ELSE
             write (iout,2) (ii,ii=ibeg,iend)
         END IF
         IF ( PRESENT(rowlab) ) THEN
              IF (frmt == 'f') then
                  do j=1,n
                     write (iout,3) rowlab(j),  (a(j,k),k=ibeg,iend)
                 end do
              else
                  do j=1,n
	             write (iout,4) rowlab(j),  (a(j,k),k=ibeg,iend)
                  end do
              end if
         ELSE
              IF (frmt == 'f') then
                  do j=1,n
                     write (iout,5) (a(j,k),k=ibeg,iend)
                  end do
              else
                  do j=1,n
	             write (iout,6)  (a(j,k),k=ibeg,iend)
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
             write(iout,1) collab(ibeg:iend)
         ELSE
             write (iout,2) (ii,ii=ibeg,iend)
         END IF
         IF ( PRESENT(rowlab) ) THEN
              do j=1,n
                 write (iout,4) rowlab(j),  (a(j,k),k=ibeg,iend)
              end do
         ELSE
              do j=1,n
                 write (iout,6) (a(j,k),k=ibeg,iend)
              end do
         END IF
         ibeg=iend
      end do
  END IF
1 format(15x,5(3x,a8,9x))
2 format(1x,'col',5(5x,i6,9x))
3 format(4x,a6,5f20.12)
4 format(4x,a6,5e20.12)
5 format(4x,5f20.12)
6 format(4x,5e20.12)
end subroutine  Print_Matrix_d
!***********************************************************************
!***********************************************************************
!deck Print_Matrix_z
  subroutine Print_Matrix_z(type_complex_matrix,a,n,m,frmt,title,collab,rowlab)
  implicit none
  TYPE(Complex_Matrix)                     :: type_complex_matrix
  COMPLEX(idp), DIMENSION(:,:)             :: a
  INTEGER                                  :: n
  INTEGER                                  :: m
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: rowlab
  INTEGER                                  :: lenth
  INTEGER                                  :: ibeg
  INTEGER                                  :: iend
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: k
  INTEGER                                  :: ii
  INTEGER                                  :: len
!
  IF (PRESENT(title)) THEN
      len=lenth(title)
      write(iout,*) title(1:len)
  END IF
  ibeg=0
  IF (PRESENT(frmt)) THEN
      do i=1,m,2
         iend=min(ibeg+2,m)
         ibeg=ibeg+1
         IF (PRESENT(collab)) THEN
             write(iout,1) collab(ibeg:iend)
         ELSE
             write (iout,2) (ii,ii=ibeg,iend)
         END IF
         IF ( PRESENT(rowlab) ) THEN
              IF (frmt == 'f') then
                  do j=1,n
                     write (iout,3) rowlab(j),  (a(j,k),k=ibeg,iend)
                 end do
              else
                  do j=1,n
	             write (iout,4) rowlab(j),  (a(j,k),k=ibeg,iend)
                  end do
              end if
         ELSE
              IF (frmt == 'f') then
                  do j=1,n
                     write (iout,5) (a(j,k),k=ibeg,iend)
                  end do
              else
                  do j=1,n
	             write (iout,6)  (a(j,k),k=ibeg,iend)
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
             write(iout,1) collab(ibeg:iend)
         ELSE
             write (iout,2) (ii,ii=ibeg,iend)
         END IF
         IF ( PRESENT(rowlab) ) THEN
              do j=1,n
                 write (iout,5) rowlab(j),  (a(j,k),k=ibeg,iend)
              end do
         ELSE
              do j=1,n
                 write (iout,6) (a(j,k),k=ibeg,iend)
              end do
         END IF
         ibeg=iend
      end do
  END IF
1 format(15x,2(15x,a8,18x))
2 format(1x,'col(Re,Im)',4(9x,i6,25x))
3 format(4x,a6,4f20.12)
4 format(4x,a6,4e20.12)
5 format(4x,4f20.12)
6 format(4x,4e20.12)
end subroutine  Print_Matrix_z
!***********************************************************************
!***********************************************************************
!deck Print_Triangle_Matrix_d
  subroutine Print_Triangle_Matrix_d(type_real_triangle,a,n,frmt,title,collab)
  implicit none
  TYPE(Real_Triangle)                      :: type_real_triangle
  REAL(idp), DIMENSION (:)                 :: a
  INTEGER                                  :: n
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  INTEGER                                  :: lenth
  INTEGER                                  :: ibeg
  INTEGER                                  :: iend
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: ii
  INTEGER                                  :: count
  INTEGER                                  :: len
!
  IF (PRESENT(title)) THEN
      len=lenth(title)
      write(iout,*) title(1:len)
  END IF
  ibeg = 0
  count = 0
  IF (PRESENT(frmt)) THEN
      DO i = 1, n, 5
         iend = min(ibeg+5,n)
         ibeg = ibeg + 1
         IF (PRESENT(collab)) THEN
             write(iout,1) collab(ibeg:iend)
         ELSE
             write (iout,2) (ii,ii=ibeg,iend)
         END IF
         IF (frmt == 'f') THEN
             DO j = ibeg, iend
                write (iout,3) a(count + 1:count + j)
                count=count+j
             END DO
         ELSE
             DO j = ibeg, iend
                write (iout,4) a(count + 1:count + j)
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
             write(iout,1) collab(ibeg:iend)
         ELSE
             write (iout,2) (ii,ii=ibeg,iend)
         END IF
         DO j = ibeg, iend
            write (iout,4) a(count + 1:count + j)
            count=count+j
         END DO
         ibeg = iend
      END DO
  END IF
1 format(15x,5(3x,a8,9x))
2 format(1x,'col',5(5x,i6,9x))
3 format(10x,5f20.12)
4 format(10x,5e20.12)
end subroutine  Print_Triangle_Matrix_d
!***********************************************************************
!***********************************************************************
!deck Print_Triangle_Matrix_z
  subroutine Print_Triangle_Matrix_z(type_complex_triangle,a,n,frmt,title,collab)
  implicit none
  TYPE(Complex_Triangle)                   :: type_complex_triangle
  COMPLEX(idp), DIMENSION (:)              :: a
  INTEGER                                  :: n
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  INTEGER                                  :: lenth
  INTEGER                                  :: ibeg
  INTEGER                                  :: iend
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: ii
  INTEGER                                  :: count
  INTEGER                                  :: len
!
  IF (PRESENT(title)) THEN
      len=lenth(title)
      write(iout,*) title(1:len)
  END IF
  ibeg = 0
  count = 0
  IF (PRESENT(frmt)) THEN
      DO i = 1, n, 5
         iend = min(ibeg+5,n)
         ibeg = ibeg + 1
         IF (PRESENT(collab)) THEN
             write(iout,1) collab(ibeg:iend)
         ELSE
             write (iout,2) (ii,ii=ibeg,iend)
         END IF
         IF (frmt == 'f') THEN
             DO j = ibeg, iend
                write (iout,3) a(count + 1:count + j)
                count=count+j
             END DO
         ELSE
             DO j = ibeg, iend
                write (iout,4) a(count + 1:count + j)
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
             write(iout,1) collab(ibeg:iend)
         ELSE
             write (iout,2) (ii,ii=ibeg,iend)
         END IF
         DO j = ibeg, iend
            write (iout,4) a(count + 1:count + j)
            count=count+j
         END DO
         ibeg = iend
      END DO
  END IF
1 format(15x,3(3x,a8,9x))
2 format(1x,'col',3(5x,i6,9x))
3 format(10X,5f20.12)
4 format(10x,5e20.12)
end subroutine  Print_Triangle_Matrix_z
!***********************************************************************
!***********************************************************************
!deck Print_Vector_d
  subroutine Print_Vector_d(type_real_vector,a,frmt,title,collab)
  implicit none
  TYPE(Real_Vector)                        :: type_real_vector
  REAL(idp), DIMENSION (:)                 :: a
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  INTEGER                                  :: lenth
  INTEGER                                  :: len
  IF (PRESENT(title)) THEN
      len=lenth(title)
      write(iout,*) title(1:len)
  END IF
  IF (PRESENT(frmt)) THEN
      IF (PRESENT(collab)) THEN
          write(iout,1) collab
      ELSE
          write (iout,2)
      END IF
      IF (frmt == 'f') THEN
          write (iout,3) a(:)
      ELSE
          write (iout,4) a(:)
      END IF
  ELSE
      IF (PRESENT(collab)) THEN
          write(iout,1) collab
      ELSE
          write (iout,2)
      END IF
      write (iout,4) a(:)
  END IF
1 format(18x,a8)
2 format(18x,'col',5x,i6)
3 format(10x,f20.12)
4 format(10x,e20.12)
end subroutine  Print_Vector_d
!***********************************************************************
!***********************************************************************
!deck Print_Vector_z
  subroutine Print_Vector_z(type_complex_vector,a,frmt,title,collab)
  implicit none
  TYPE(Complex_Vector)                     :: type_complex_vector
  COMPLEX(idp), DIMENSION (:)              :: a
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  INTEGER                                  :: lenth
  INTEGER                                  :: len
!
  IF (PRESENT(title)) THEN
      len=lenth(title)
      write(iout,*) title(1:len)
  END IF
  IF (PRESENT(frmt)) THEN
      IF (PRESENT(collab)) THEN
          write(iout,1) collab
      ELSE
          write (iout,2)
      END IF
      IF (frmt == 'f') THEN
          write (iout,3) a(:)
      ELSE
          write (iout,4) a(:)
      END IF
  ELSE
      IF (PRESENT(collab)) THEN
          write(iout,1) collab
      ELSE
          write (iout,2)
      END IF
      write (iout,4) a(:)
  END IF
1 format(18x,a8)
2 format(18x,'col',5x,i6)
3 format(10X,2f20.12)
4 format(10x,2e20.12)
end subroutine  Print_Vector_z
!***********************************************************************
!***********************************************************************
END MODULE Matrix_Print
