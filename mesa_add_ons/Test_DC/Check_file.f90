
!======================================================================
      Subroutine Check_file(AF)
!======================================================================
      Character(*), Intent(in) :: AF
      Logical :: EX
      Inquire(file=AF,exist=EX)
      if(.not.EX) then
       write(*,*) ' can not find file  ',AF
       Stop
      end if
      End Subroutine Check_file
