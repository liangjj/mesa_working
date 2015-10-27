!======================================================================
      Subroutine Read_apar(nu,name,avalue)
!======================================================================
!
!     read the character variable 'avalue' with identifier 'name'
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
      Character(*) :: avalue
      Character(80) :: AS
      Integer(4) :: i
      i=LEN_TRIM(name)
      rewind(nu)
    1 read(nu,'(a)',end=2) AS
    if(AS(1:i).ne.name) go to 1
      i=INDEX(AS,'=')+1
      read(AS(i:),*) avalue
    2 Continue
      End Subroutine Read_apar
