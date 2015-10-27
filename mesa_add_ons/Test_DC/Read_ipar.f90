

!======================================================================
      Subroutine Read_ipar(nu,name,ivalue)
!======================================================================
!
!     read the integer variable 'ivalue' with identifier 'name'
!     from unit 'nu', where the record like
!
!     name =  #####      is supposed to exist
!
!----------------------------------------------------------------------
      Implicit none
      Integer(4), Intent(in) :: nu
      Character(*), Intent(in) :: name
      Integer(4) :: ivalue
      Character(80) :: AS
      Integer(4) :: i
      i=LEN_TRIM(name)
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
    if(AS(1:i).ne.name) go to 1
      i=INDEX(AS,'=')+1
      read(AS(i:),*) ivalue
    2 Continue
      End Subroutine Read_ipar
