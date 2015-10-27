PROGRAM test
   IMPLICIT NONE
   
   INTEGER, PARAMETER :: UnitOut = 6
   INTEGER, ALLOCATABLE, DIMENSION(:) :: GroupListCnt
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: GroupListMem

   CHARACTER (LEN = 5), ALLOCATABLE, DIMENSION(:,:) :: GroupListSym

   INTEGER :: nsmall
   INTEGER :: NGroupOcc
   INTEGER :: i, j

   nsmall = 11
   NGroupOcc = 2
   

   ALLOCATE (GroupListMem(nsmall, ABS(NGroupOcc)))
   ALLOCATE (GroupListSym(nsmall, ABS(NGroupOcc)))
   ALLOCATE (GroupListCnt(ABS(NGroupOcc)))
   
   GroupListCnt(1) = 10
   GroupListCnt(2) = 1

   GroupListMem = 1
   GroupListSym( 1, 1) = 'a    '
   GroupListSym( 2, 1) = 'b    '
   GroupListSym( 3, 1) = 'c    '
   GroupListSym( 4, 1) = 'd    '
   GroupListSym( 5, 1) = 'e    '
   GroupListSym( 6, 1) = 'f    '
   GroupListSym( 7, 1) = 'g    '
   GroupListSym( 8, 1) = 'h    '
   GroupListSym( 9, 1) = 'i    '
   GroupListSym(10, 1) = 'j    '

   GroupListSym( 1, 2) = 'k    '

   DO i = 1, ABS(NGroupOcc)
      WRITE (UNIT = UnitOut, FMT = "('Group number ', i5, '  with', i5, '  orbitals')") i, GroupListCnt(i)
      WRITE (UNIT = UnitOut, FMT = "((5(i5, '-', a5)))") (GroupListMem(j, i), GroupListSym(j, i),&
              &  j = 1, GroupListCnt(i))
   END DO
   DO i = 1, ABS(NGroupOcc)
      WRITE (UNIT = UnitOut, FMT = "('Group number ', i5, '  with', i5, '  orbitals')") i, GroupListCnt(i)
      WRITE (UNIT = UnitOut, FMT = "((5(i5, '-', a5)))") (GroupListMem(j, i), ADJUSTL(GroupListSym(j, i)),&
              &  j = 1, GroupListCnt(i))
   END DO
END PROGRAM test
