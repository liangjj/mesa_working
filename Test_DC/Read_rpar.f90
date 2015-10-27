!======================================================================
      Subroutine Read_rpar(nu,name,rvalue)
!======================================================================
!
!     read the real variable 'rvalue' with identifier 'name'
!     from unit 'nu', where the record like
!
!     name =  #####    ! coments
!
!     is supposed to exist
!
!----------------------------------------------------------------------
      Implicit none
      Integer(4), Intent(in) :: nu
      Character(*), Intent(in) :: name
      Real(8) :: rvalue
      Character(80) :: AS
      Integer(4) :: i
      i=LEN_TRIM(name)
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
    if(AS(1:i).ne.name) go to 1
      i=INDEX(AS,'=')+1
      read(AS(i:),*) rvalue
    2 Continue
      End Subroutine Read_rpar
