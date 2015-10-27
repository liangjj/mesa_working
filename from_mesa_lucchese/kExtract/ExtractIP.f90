! ** 67890123 20 890123 30 890123 40 890123 50 890123 60 890123 70 890123 80 890123 90 890123 100 90123 110 90123 120 90
PROGRAM ExtractIP
! ** rad configuration lists from mesa and count the number of such terms
   USE Extract_mod

   IMPLICIT NONE

   INTEGER (KIND = 4), PARAMETER :: UnitOutIP = 10, UnitIn = 11, UnitOut = 12, UnitSh = 13, UnitTmp = 14
   CHARACTER (LEN = InpLineL) :: Line
   INTEGER (KIND = 4) :: ios

   INTEGER :: NChan
   REAL (KIND = XR) :: FirstIP
   REAL (KIND = XR), ALLOCATABLE, DIMENSION(:) :: IPot
   INTEGER :: i

   OPEN (UNIT = UnitIn, FILE = "inExtractIP", POSITION = 'REWIND')
   OPEN (UNIT = UnitOut, FILE = "outExtractIP", POSITION = 'REWIND')
   
   WRITE (UNIT = UnitOut, FMT ="('Program ExtractIP - get ionization potentials')")

   REWIND (UNIT = UnitIn)
   READ (UNIT = UnitIn, FMT = *) NChan
   READ (UNIT = UnitIn, FMT = *) FirstIP

   WRITE (UNIT = UnitOut, FMT = "('Number of channels' , i5)") NChan
   WRITE (UNIT = UnitOut, FMT = "('Ionization potential in the first channel', e17.8)") FirstIP
   ALLOCATE (IPot(NChan))

   OPEN (UNIT = UnitOutIP, FILE = "mesa.out", POSITION = 'REWIND')

   DO
      READ (UNIT = UnitOutIP, IOSTAT = ios, FMT = "(a)") Line
      IF (ios /= 0) THEN
         STOP 'Could not find Input Eigenvalues'
      END IF
      IF (Line == "  Input Eigenvalues")  THEN
         EXIT
      END IF
   END DO

   READ (UNIT = UnitOutIP, FMT = *) IPot
   WRITE (UNIT = UnitOut, FMT = "('Input eigenvalues', /, (5f15.6))") IPot
   IPot = (IPot-IPot(1))*XFEVAU + FirstIP
   WRITE (UNIT = UnitOut, FMT = "('Ionization potentials for inphot', /, (5f15.6))") IPot

   OPEN (UNIT = UnitTmp, FILE = "inphot", POSITION = 'REWIND')
   WRITE (UNIT = UnitTmp, FMT = "(3I5)") 0, 0   ! iprt, nww - write out all three transition moments
   WRITE (UNIT = UnitTmp, FMT = "(10f12.6)") IPot
   CLOSE (UNIT = UnitTmp)

   OPEN (UNIT = UnitSh, FILE = "ExtractIP.sh", POSITION = 'REWIND')
   WRITE (UNIT = UnitSh, FMT = "('#')")
   WRITE (UNIT = UnitSh, FMT = "('IPot=( ')")
   DO i = 1, NChan
      WRITE (UNIT = UnitSh, FMT = "(f12.6)") IPot(i)
   END DO
   WRITE (UNIT = UnitSh, FMT = "(')')")
   
   CLOSE (UNIT = UnitSh)

END PROGRAM ExtractIP
