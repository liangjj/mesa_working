PROGRAM ExtractData
   USE Extract_mod

   IMPLICIT NONE

   INTEGER (KIND = 4), PARAMETER :: UnitOutScf = 10, UnitIn = 11, UnitOut = 12, UnitSh = 13, UnitTmp = 14, UnitTmp2 = 15
   INTEGER (KIND = 4), PARAMETER :: UnitTmp3 = 16
   INTEGER, PARAMETER :: LineCntX = 7
   INTEGER, PARAMETER :: NSymMax = 40

   REAL (KIND = XR), ALLOCATABLE, DIMENSION(:,:) :: Coords
   REAL (KIND = XR), DIMENSION(3) :: Origin
   INTEGER :: NCenter
   INTEGER :: CenterExpand 

   CHARACTER (LEN = InpLineL) :: Line
   INTEGER :: LineCnt
   INTEGER (KIND = 4) :: ios
   INTEGER :: icen, jcen, i, L, M, j, k, ismall, iscat
   INTEGER :: NSym
   CHARACTER (LEN = 5), ALLOCATABLE, DIMENSION(:) :: SymLabel
   CHARACTER (LEN = 5), ALLOCATABLE, DIMENSION(:) :: SymLabelePoly
   CHARACTER (LEN = 5) :: SymIn
   INTEGER, ALLOCATABLE, DIMENSION(:) :: SymNumOrb
   INTEGER, PARAMETER :: NBuff = 500

   TYPE (IntVec), POINTER, DIMENSION(:) :: FFBFData, KOHNOPTData
   TYPE (IntVec), POINTER, DIMENSION(:) :: OrbInSmall, OrbInScat
   INTEGER, ALLOCATABLE, DIMENSION(:) :: NOrbInSmall, NOrbInScat
   INTEGER, ALLOCATABLE, DIMENSION(:) :: FFBFN
   INTEGER :: KOHNOPTN
   INTEGER :: nsmall
   CHARACTER (LEN = 10), ALLOCATABLE, DIMENSION(:) :: SymCont  ! symmetry label for the continuum 1-electron functions
   CHARACTER (LEN = 10) :: SymScat  ! symmetry label of the total N-electron scattering state
   CHARACTER (LEN = 10) :: SymInit  ! symmetry label of the total N-electron initial state state
   CHARACTER (LEN = 10), ALLOCATABLE, DIMENSION(:) :: SymTarg  ! symmetry label of the total N-1 electron ion state
   INTEGER, ALLOCATABLE, DIMENSION(:) :: nscont, nstarg! numeric representation for the Cont and Scat symmetries
   INTEGER :: nsscat, nsinit
   INTEGER :: LMaxK, CenterSym, CnOrder
   INTEGER, ALLOCATABLE, DIMENSION(:) :: PWLStart, PWmMOD, iCS
   INTEGER, ALLOCATABLE, DIMENSION(:) :: npw
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: PWList
   INTEGER :: ScatGroup
   INTEGER :: Symmetry ! number of different symmetries in the the scattering channels <= nchan
   INTEGER :: iSym   
   INTEGER :: nchan   ! number of scattering channels
   INTEGER :: ichan   ! index for nchan
   INTEGER, ALLOCATABLE, DIMENSION(:) :: symroot
   INTEGER :: NGroupOcc
   INTEGER, ALLOCATABLE, DIMENSION(:) :: GroupListCnt
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: GroupListMem
   CHARACTER (LEN = 5), ALLOCATABLE, DIMENSION(:,:) :: GroupListSym
   INTEGER, ALLOCATABLE, DIMENSION(:) :: GroupList
   INTEGER :: na, nb ! defines the spin state for mesa
   INTEGER, ALLOCATABLE, DIMENSION(:) :: OrbSymList
   CHARACTER (LEN = 9) :: iconq

   OPEN (UNIT = UnitIn, FILE = "inExtractData", POSITION = 'REWIND')
   OPEN (UNIT = UnitOut, FILE = "outExtractData", POSITION = 'REWIND')

   WRITE (UNIT = UnitOut, FMT = "('Program ExtractData - get data from mesa.out')")

! ** read input from inExtractData

   REWIND (UNIT = UnitIn)

   READ (UNIT = UnitIn, FMT = *) nsmall
   WRITE (UNIT = UnitOut, FMT = "('nsmall = ', i5)") nsmall
   READ (UNIT = UnitIn, FMT = *) na, nb
   WRITE (UNIT = UnitOut, FMT = "('na', i5, '  nb', i5)") na, nb
   READ (UNIT = UnitIn, FMT = *) CenterExpand
   WRITE (UNIT = UnitOut, FMT = "('Atomic center to use for the asymptotic functions ', i5)") CenterExpand
                                ! if zero then leave expansion center where it is
   READ (UNIT = UnitIn, FMT = *) iconq
   WRITE (UNIT = UnitOut, FMT = "('Label for insolve for contractq is ', a)") iconq
   READ (UNIT = UnitIn, FMT = *) LMaxK
   WRITE (UNIT = UnitOut, FMT = "('Highest asymptotic partial wave LMaxK ', i5)") LMaxK
   READ (UNIT = UnitIn, FMT = *) CenterSym, CnOrder  
                                ! CenterSym = 1 for no center of symmetry, = 2 with center of symmetry
                                ! CnOrder = is the highest rotation symmetry about the z axis
   WRITE (UNIT = UnitOut, FMT = "('Center of symmeter (1 no, 2 yes) CenterSym ', i5)") CenterSym
   WRITE (UNIT = UnitOut, FMT = "('Rotation symmetry order about z axis CnOrder', i5)") CnOrder

   READ (UNIT = UnitIn, FMT = *) nchan, Symmetry
   WRITE (UNIT = UnitOut, FMT = "('Number of scattering channels', i5)") nchan
   WRITE (UNIT = UnitOut, FMT = "('Number of continuum symmetry types' , i5)") Symmetry
   ALLOCATE (symroot(nchan))
   READ (UNIT = UnitIn, FMT = *) symroot
   WRITE (UNIT = UnitOut, FMT = "('Symmetry type of each channel', 10i5)") symroot
   READ (UNIT = UnitIn, FMT = *) NGroupOcc  !  number of groups used for the occupied orbitals in the DRTs
                                ! If NGroupOcc is positive then read in orbitals number in the total list
                                ! If NGroupOcc is negative then read in both orbital number and symmetry type
                                ! where the orbital numbers are then for each symmetry type, e.g. 4 'a1'
                                ! for the 4th a1 orbital
   WRITE (UNIT = UnitOut, FMT = "('Number of diffrerent groups for nsmall orbitals NGroupOcc', i5)") NGroupOcc
   IF (NGroupOcc == 0) STOP "Bad nGroupOcc"
   ALLOCATE (GroupListMem(nsmall, ABS(NGroupOcc)))
   ALLOCATE (GroupListCnt(ABS(NGroupOcc)))
   IF (NGroupOcc > 0) THEN
      DO i = 1, NGroupOcc
         READ (UNIT = UnitIn, FMT = *) GroupListCnt(i), (GroupListMem(j, i), j = 1, GroupListCnt(i))
      END DO
   ELSE
      ALLOCATE (GroupListSym(nsmall, ABS(NGroupOcc)))
      DO i = 1, ABS(NGroupOcc)
         READ (UNIT = UnitIn, FMT = *) GroupListCnt(i), (GroupListMem(j, i), GroupListSym(j, i), j = 1, GroupListCnt(i))
         IF (GroupListCnt(i) > nsmall) STOP 'nsmall not large enough'
         DO j = 1, GroupListCnt(i)
            GroupListSym(j, i) = ADJUSTL(GroupListSym(j, i))
         END DO
      END DO
   END IF
   DO i = 1, ABS(NGroupOcc)
      WRITE (UNIT = UnitOut, FMT = "('Group number ', i5, '  with', i5, '  orbitals')") i, GroupListCnt(i)
      IF (NGroupOcc > 0) THEN
         WRITE (UNIT = UnitOut, FMT = "(10i5)") (GroupListMem(j, i), j = 1, GroupListCnt(i))
      ELSE
         WRITE (UNIT = UnitOut, FMT = "((5(i5, '-', a5)))") (GroupListMem(j, i), GroupListSym(j, i),&
              &  j = 1, GroupListCnt(i))
      END IF
   END DO


   ALLOCATE (SymCont(Symmetry), SymTarg(Symmetry))
   ALLOCATE (PWLStart(Symmetry), PWmMOD(Symmetry), iCS(Symmetry))
   READ (UNIT = UnitIn, FMT = *) SymInit
   WRITE (UNIT = UnitOut, FMT = "('Symmetry of the total N electron initial state ', a)") &
        & TRIM(ADJUSTL(SymInit))
   READ (UNIT = UnitIn, FMT = *) SymScat
   WRITE (UNIT = UnitOut, FMT = "('Symmetry of the total N electron scattering state ', a)") &
        & TRIM(ADJUSTL(SymScat))
   DO iSym = 1, Symmetry
      WRITE (UNIT = UnitOut, FMT = "('For Symmetry type', i5)") iSym
      READ (UNIT = UnitIn, FMT = *) SymCont(iSym), SymTarg(iSym)
      WRITE (UNIT = UnitOut, FMT = "('Symmetry for the continuum one electron functions ', a)") &
           & TRIM(ADJUSTL(SymCont(iSym)))
      WRITE (UNIT = UnitOut, FMT = "('Symmetry of the total N-1 electron+target state ', a)") &
           & TRIM(ADJUSTL(SymTarg(iSym)))
      READ (UNIT = UnitIn, FMT = *) PWLStart(iSym), PWmMOD(iSym), iCS(iSym)  ! for this set of asymptotic waves
                                !   PWLStart = 0 for even wave when there is no even-odd symmetry
                                !   PWLStart = 1 for only odd waves
                                !   PWmMOD, only include wave when MOD(m, CnOrder) = PWmMOD
                                !   iCS = 0 for cos(phi) like terms, = 1 for sin(phi) like terms
      WRITE (UNIT = UnitOut, FMT = "('Lowest L for asymptotic partial waves PWLStart', i5)") PWLStart(iSym)
      WRITE (UNIT = UnitOut, FMT = "('MOD of M values to use PWmMOD', i5)") PWmMOD(iSym)
      WRITE (UNIT = UnitOut, FMT = "('Cos (0) or Sin (1) like terms iCS', i5)") iCS(iSym)
   END DO


   ScatGroup = ABS(NGroupOcc) + 1

! ** read the output from the scf run in mesa.out

   OPEN (UNIT = UnitOutScf, FILE = "mesa.out", POSITION = 'REWIND')

   REWIND (UNIT = UnitOutScf)

   DO
      READ (UNIT = UnitOutScf, FMT = "(a)", IOSTAT = ios) Line
      IF (ios /= 0) THEN
         STOP "could not find geometry"
      END IF
      IF (Line == "     cd cent  el                   coordinates(Bohr)") THEN
         EXIT
      END IF
   END DO

   READ (UNIT = UnitOutScf, FMT = "(a)", IOSTAT = ios) Line
   IF (ios /= 0) THEN
      STOP "could not find geometry xyz"
   END IF

   NCenter = 0
   DO
      READ (UNIT = UnitOutScf, FMT = "(a)", IOSTAT = ios) Line
      IF (ios /= 0) THEN
         STOP "could not find geometry lines"
      END IF
      IF (Line(1:28) == " rotational constants (ghz):") THEN
         EXIT
      END IF
      READ (UNIT = Line, FMT = "(i7, i4)") jcen, icen
      IF (jcen /= 0) THEN
         NCenter = NCenter + 1
      END IF
   END DO

   IF (NCenter == 0) STOP "could not find centers"

   ALLOCATE (Coords(3, NCenter))

   REWIND (UNIT = UnitOutScf)

   DO
      READ (UNIT = UnitOutScf, FMT = "(a)", IOSTAT = ios) Line
      IF (ios /= 0) THEN
         STOP "could not find geometry"
      END IF
      IF (Line == "     cd cent  el                   coordinates(Bohr)") THEN
         EXIT
      END IF
   END DO

   READ (UNIT = UnitOutScf, FMT = "(a)", IOSTAT = ios) Line
   IF (ios /= 0) THEN
      STOP "could not find geometry xyz"
   END IF

   NCenter = 0
   DO
      READ (UNIT = UnitOutScf, FMT = "(a)", IOSTAT = ios) Line
      IF (ios /= 0) THEN
         STOP "could not find geometry lines"
      END IF
      IF (Line(1:28) == " rotational constants (ghz):") THEN
         EXIT
      END IF
      READ (UNIT = Line, FMT = "(i7, i4)") jcen, icen
      IF (jcen /= 0) THEN
         NCenter = NCenter + 1
         READ (UNIT = Line, FMT = "(25x,3f12.6)") Coords(:,NCenter)
! **         10        20        30        40        50        60        70
! ** 1234567890123456789012345678901234567890123456789012345678901234567890
! **       1   1   n              0.000000    0.000000    0.615463
      END IF
   END DO


   REWIND (UNIT = UnitOutScf)

   DO
      READ (UNIT = UnitOutScf, FMT = "(a)", IOSTAT = ios) Line
      IF (ios /= 0) THEN
         STOP "could not find symmetry"
      END IF
! **         10        20        30        40        50        60        70
! ** 1234567890123456789012345678901234567890123456789012345678901234567890
! **      number of salcs :   a1   a2   b1   b2   
! **                          34   10   13   25

      IF (Line(1:22) == "     number of salcs :") THEN
         EXIT
      END IF
   END DO

   DO NSym = 1, NSymMax
      READ (UNIT = Line, FMT = "(22x, 3x, 40a5)") (SymIn, i = 1, NSym)
      IF (SymIn == " ") THEN
         EXIT
      END IF
   END DO
   IF (NSym > NSymMax) STOP "NSymMax not large enough"
   NSym = NSym -1
   IF (NSym == 0) STOP "Did not find any symmetries"
   WRITE (UNIT = UnitOut, FMT = "('NSym =', i5)") NSym

   ALLOCATE (SymLabel(NSym))
   ALLOCATE (SymLabelePoly(NSym))
   ALLOCATE (SymNumOrb(NSym))
   READ (UNIT = Line, FMT = "(22x, 3x, 40a5)") SymLabel
   READ (UNIT = UnitOutScf, FMT = "(22x, 40i5)", IOSTAT = ios) SymNumOrb
   IF (ios /= 0) THEN
      STOP "could not find number of orbitals of each symmetry"
   END IF

   SymLabel = ADJUSTL(SymLabel)
   WRITE (UNIT = UnitOut, FMT = "('SymLabel   =', 10a5)") SymLabel
   WRITE (UNIT = UnitOut, FMT = "('SymNormOrb =', 10i5)") SymNumOrb

   DO i = 1, NSym
      CALL ToCaps(SymLabel(i), SymLabelePoly(i))
   END DO

   REWIND (UNIT = UnitOutScf)
   DO
      READ (UNIT = UnitOutScf, FMT = "(a)", IOSTAT = ios) Line
      IF (ios /= 0) THEN
         STOP "could not FFBF Input"
      END IF
! **         10        20        30        40        50        60        70
! ** 1234567890123456789012345678901234567890123456789012345678901234567890
! **      number of salcs :   a1   a2   b1   b2   
! **                          34   10   13   25

      IF (Line(1:22) == " FFBF Input") THEN
         EXIT
      END IF
   END DO

   ALLOCATE (FFBFData(NSym), FFBFN(NSym), KOHNOPTData(NSym))

   DO i = 1, NSym
      READ (UNIT = UnitOutScf, FMT = "()")
      READ (UNIT = UnitOutScf, FMT = "()")
      READ (UNIT = UnitOutScf, FMT = *) FFBFN(i)
      ALLOCATE (FFBFData(i)%v(FFBFN(i)))
      READ (UNIT = UnitOutScf, FMT = *) FFBFData(i)%v
   END DO

   READ (UNIT = UnitOutScf, FMT = "()")
   READ (UNIT = UnitOutScf, FMT = "(a)") Line
   IF (Line /= " KOHNOPT Input") STOP "KOHNOPT not where it is expected"

   DO i = 1, NSym
      READ (UNIT = UnitOutScf, FMT = "()")
      READ (UNIT = UnitOutScf, FMT = "()")
      READ (UNIT = UnitOutScf, FMT = *) KOHNOPTN
      IF (KOHNOPTN /= SymNumOrb(i)) STOP "KOHNOPTN not the expected size"
      ALLOCATE (KOHNOPTData(i)%v(KOHNOPTN))
      READ (UNIT = UnitOutScf, FMT = *) KOHNOPTData(i)%v
   END DO
   CLOSE (UNIT = UnitOutScf)

   WRITE (UNIT = UnitOut, FMT = "('analyze OrbInSmall NSym =', i5)") NSym
! ** divide the molecular orbital lists for those in nsmall and those not in nsmall

   ALLOCATE (OrbInSmall(NSym), NOrbInSmall(NSym), OrbInScat(NSym), NOrbInScat(NSym))
   ALLOCATE (OrbSymList(nsmall))
   OrbSymList = 0
   NOrbInSmall = 0
   NOrbInScat = 0
   DO i = 1, NSym
      NOrbInSmall(i)= COUNT(KOHNOPTData(i)%v <= nsmall) 
      NOrbInScat(i)= COUNT(KOHNOPTData(i)%v > nsmall)
      WRITE (UNIT = UnitOut, FMT = "('iSym =', i5, '  NOrbInSmall =', i5, '  NOrbInScat =', i5)") i, NOrbInSmall(i),&
           & NOrbInScat(i)

      ALLOCATE (OrbInSmall(i)%v(NOrbInSmall(i)))
      ALLOCATE (OrbInScat(i)%v(NOrbInScat(i)))
      ismall = 0
      iscat = 0
      DO j = 1, SymNumOrb(i)
         IF (KOHNOPTData(i)%v(j) <= nsmall) THEN
            ismall = ismall + 1
            OrbInSmall(i)%v(ismall) = KOHNOPTData(i)%v(j)
            OrbSymList(KOHNOPTData(i)%v(j)) = i
         ELSE
            iscat = iscat + 1
            OrbInScat(i)%v(iscat) = KOHNOPTData(i)%v(j)
         END IF
      END DO
   END DO
   IF (ANY(OrbSymList == 0)) STOP "Symmetry of orbital not identified"

   WRITE (UNIT = UnitOut, FMT = "(/, 'Orbitals numbers by symmetry', /, '  Sym Num    Orb Num')")
   DO i = 1, NSym
      DO j = 1, NOrbInSmall(i)
         WRITE (UNIT = UnitOut, FMT = "(i5, '-', a5, i7)") j, SymLabel(i), OrbInSmall(i)%v(j)
      END DO
   END DO

   WRITE (UNIT = UnitOut, FMT = "(/, 'Orbitals by number', /, 'Orb Num  Sym Num')")
   DO k = 1, NSmall
      i = OrbSymList(k)
      DO j = 1, NOrbInSmall(i)
         IF (k == OrbInSmall(i)%v(j)) THEN
            EXIT
         END IF
      END DO
      IF (j > NOrbInSmall(i)) STOP  'Inconsistent Orbital Lists'
      WRITE (UNIT = UnitOut, FMT = "(i5, i7, '-', a5)") k, j, SymLabel(i)
   END DO
   WRITE (UNIT = UnitOut, FMT = "()")

! ** Determine symmetry type for the scattering state and continuum functions

   nsscat = GetSymTyp(SymScat, SymLabel)
   nsinit = GetSymTyp(SymInit, SymLabel)
   
   ALLOCATE (nscont(Symmetry), nstarg(Symmetry))
   DO iSym = 1, Symmetry
      nscont(iSym) = GetSymTyp(SymCont(iSym), SymLabel)
      nstarg(iSym) = GetSymTyp(SymTarg(iSym), SymLabel)
   END DO

! ** get actual orbital numbers for DRT groups
! ** GroupList(i) is the group number for the i'th orbital
   
   ALLOCATE (GroupList(nsmall))
   GroupList = 0
   DO i = 1, ABS(NGroupOcc)
      DO j = 1, GroupListCnt(i)
         IF (NGroupOcc > 0) THEN
            GroupList(GroupListMem(j, i)) = i
         ELSE
            iSym = GetSymTyp(GroupListSym(j, i), SymLabel)
            IF (GroupListMem(j, i) > NOrbInSmall(iSym)) STOP "Bad GroupListMem value"
            GroupList(OrbInSmall(iSym)%v(GroupListMem(j, i))) = i
         END IF
      END DO
   END DO
   IF (ANY(GroupList == 0)) STOP "Orbital not in group list"

! ** drt for set

   OPEN (UNIT = UnitTmp, FORM = "FORMATTED", POSITION = "REWIND", FILE = "setdrt")
   WRITE (UNIT = UnitTmp, FMT = "('$drt')")
   CALL WriteOccDrt
   DO iSym = 1, Symmetry
      CALL AddToDRTLine(Line, 1, ScatGroup+iSym-1, nscont(iSym), LineCnt, LineCntX, UnitTmp)
   END DO

   CALL AddToDRTLine(Line, SUM(NOrbInScat)-Symmetry, ScatGroup+Symmetry, 1, LineCnt, LineCntX, UnitTmp)
   IF (LineCnt > 0) THEN
      WRITE (UNIT = UnitTmp, FMT = "(a)") TRIM(Line)
   END IF
   WRITE (UNIT = UnitTmp, FMT = "(a)") "na="//TRIM(Int4ToCharL(na))//" nb="//TRIM(Int4ToCharL(nb))&
        &//" ns="//TRIM(Int4ToCharL(nsscat))   
   WRITE (UNIT = UnitTmp, FMT = "('$end')")
   CLOSE (UNIT = UnitTmp)

! ** sdrt for set

   OPEN (UNIT = UnitTmp, FORM = "FORMATTED", POSITION = "REWIND", FILE = "setsdrt")
   WRITE (UNIT = UnitTmp, FMT = "('$sdrt')")
   CALL WriteOccDrt

   CALL AddToDRTLine(Line, 1, ScatGroup, 1, LineCnt, LineCntX, UnitTmp)

   CALL AddToDRTLine(Line, SUM(NOrbInScat)-1, ScatGroup+1, 1, LineCnt, LineCntX, UnitTmp)
   IF (LineCNt > 0) THEN
      WRITE (UNIT = UnitTmp, FMT = "(a)") TRIM(Line)
   END IF
   WRITE (UNIT = UnitTmp, FMT = "(a)") "na="//TRIM(Int4ToCharL(na))//" nb="//TRIM(Int4ToCharL(nb))&
        &//" ns="//TRIM(Int4ToCharL(nstarg(1)))   
   WRITE (UNIT = UnitTmp, FMT = "('$end')")
   CLOSE (UNIT = UnitTmp)

! ** sdrt for phot

   OPEN (UNIT = UnitTmp, FORM = "FORMATTED", POSITION = "REWIND", FILE = "photsdrt")
   WRITE (UNIT = UnitTmp, FMT = "('$sdrt')")
   CALL WriteOccDrt
   DO iSym = 1, Symmetry
      CALL AddToDRTLine(Line, 1, ScatGroup+iSym-1, nscont(iSym), LineCnt, LineCntX, UnitTmp)
   END DO

   CALL AddToDRTLine(Line, SUM(NOrbInScat)-Symmetry, ScatGroup+Symmetry, 1, LineCnt, LineCntX, UnitTmp)
   IF (LineCNt > 0) THEN
      WRITE (UNIT = UnitTmp, FMT = "(a)") TRIM(Line)
   END IF
   WRITE (UNIT = UnitTmp, FMT = "(a)") "na="//TRIM(Int4ToCharL(na))//" nb="//TRIM(Int4ToCharL(nb))&
        &//" ns="//TRIM(Int4ToCharL(nsscat))   
   WRITE (UNIT = UnitTmp, FMT = "('$end')")
   CLOSE (UNIT = UnitTmp)

! ** drt for opt

   OPEN (UNIT = UnitTmp, FORM = "FORMATTED", POSITION = "REWIND", FILE = "optsdrt")
   WRITE (UNIT = UnitTmp, FMT = "('$sdrt')")
   CALL WriteOccDrt

   j = 0
   DO i = 1, NSym
      IF (ALL(i /= nscont)) THEN
         j = j+NOrbInScat(i)
      ELSE
         DO k = 1, Symmetry
            IF (i == nscont(k)) THEN
               EXIT
            END IF
         END DO
         IF (j /= 0) THEN
            CALL AddToDRTLine(Line, j, ScatGroup+Symmetry, 1, LineCnt, LineCntX, UnitTmp)
            j = 0
         END IF
         CALL AddToDRTLine(Line, NOrbInScat(i), ScatGroup+k-1, i, LineCnt, LineCntX, UnitTmp)
      END IF
   END DO
   IF (j > 0) THEN
      CALL AddToDRTLine(Line, j, ScatGroup+Symmetry, 1, LineCnt, LineCntX, UnitTmp)
   END IF

   IF (LineCnt > 0) THEN
      WRITE (UNIT = UnitTmp, FMT = "(a)") TRIM(Line)
      LineCnt = 0
      Line = " "
   END IF

   WRITE (UNIT = UnitTmp, FMT = "(a)") "na="//TRIM(Int4ToCharL(na))//" nb="//TRIM(Int4ToCharL(nb))&
        &//" ns="//TRIM(Int4ToCharL(nsscat))   
   WRITE (UNIT = UnitTmp, FMT = "('$end')")
   CLOSE (UNIT = UnitTmp)

! ** write out environment variables for the mesa input

   OPEN (UNIT = UnitSh, FILE = "ExtractData.sh", POSITION = 'REWIND')
   WRITE (UNIT = UnitSh, FMT = "('#')")
   WRITE (UNIT = UnitSh, FMT = "('nsscat=', a , '  # symmetry type number for the overall scattering state')")&
        & TRIM(Int4ToCharL(nsscat))
   WRITE (UNIT = UnitSh, FMT = "('SymScat=', a)") TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "('symnl2=""', a , '""  # number for the continuum symmetry basis functions')")&
        & TRIM(MakeList(NOrbInScat(nscont)))
   WRITE (UNIT = UnitSh, FMT = "('symroot=""', a, '""  # symmetry types of the channels')") &
        & TRIM(MakeList(symroot))

   WRITE (UNIT = UnitSh, FMT = "('ngroupss=', a, '  # ngroups with different symmetry continua')") &
        & TRIM(Int4ToCharL(ScatGroup+Symmetry))
   WRITE (UNIT = UnitSh, FMT = "('ngroupsn=', a, '  # ngroups with only one continuum type')") &
        & TRIM(Int4ToCharL(ScatGroup+1))

   WRITE (UNIT = UnitSh, FMT = "('# scripts for Kohn part of the calculation')")
   WRITE (UNIT = UnitSh, FMT = "(a)") "function kohnsolve {"

   WRITE (UNIT = UnitSh, FMT = "(a)") "echo start xlebgrid ; date"
   WRITE (UNIT = UnitSh, FMT = "(a)") "xlebgrid"
   WRITE (UNIT = UnitSh, FMT = "(a)") "if [ ""$debug"" != """" ] ; then cat outgrid ; fi"
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv outgrid outgrid."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv ingrid ingrid."//TRIM(ADJUSTL(SymScat))

   WRITE (UNIT = UnitSh, FMT = "(a)") "echo start xcbetter1 ; date"
   WRITE (UNIT = UnitSh, FMT = "(a)") "xcbetter1"
   WRITE (UNIT = UnitSh, FMT = "(a)") "if [ ""$debug"" != """" ] ; then cat outcoul ; fi"
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv outcoul outcoul."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv incoul incoul."//TRIM(ADJUSTL(SymScat))

   WRITE (UNIT = UnitSh, FMT = "(a)") "echo start xass ; date"
   WRITE (UNIT = UnitSh, FMT = "(a)") "xass"
   WRITE (UNIT = UnitSh, FMT = "(a)") "if [ ""$debug"" != """" ] ; then cat outass ; fi"
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv outass outass."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv inass inass."//TRIM(ADJUSTL(SymScat))

   WRITE (UNIT = UnitSh, FMT = "(a)") "echo start xbasis ; date"
   WRITE (UNIT = UnitSh, FMT = "(a)") "xbasis"
   WRITE (UNIT = UnitSh, FMT = "(a)") "if [ ""$debug"" != """" ] ; then cat outbas ; fi"
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv outbas outbas."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv inbas inbas."//TRIM(ADJUSTL(SymScat))

   WRITE (UNIT = UnitSh, FMT = "(a)") "echo start xstatnuc ; date"
   WRITE (UNIT = UnitSh, FMT = "(a)") "xstatnuc"
   WRITE (UNIT = UnitSh, FMT = "(a)") "if [ ""$debug"" != """" ] ; then cat outstat ; fi"
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv outstat outstat."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv instat instat."//TRIM(ADJUSTL(SymScat))

   WRITE (UNIT = UnitSh, FMT = "(a)") "echo start xcffbfcg ; date"
   WRITE (UNIT = UnitSh, FMT = "(a)") "xcffbfcg"
   WRITE (UNIT = UnitSh, FMT = "(a)") "if [ ""$debug"" != """" ] ; then cat outffbf ; fi"
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv outffbf outffbf."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv inffbf inffbf."//TRIM(ADJUSTL(SymScat))

   WRITE (UNIT = UnitSh, FMT = "(a)") "echo start xcknew ; date"
   WRITE (UNIT = UnitSh, FMT = "(a)") "cat inkohnb >inkohn"
   WRITE (UNIT = UnitSh, FMT = "(a)") "xcknew"
   WRITE (UNIT = UnitSh, FMT = "(a)") "if [ ""$debug"" != """" ] ; then cat outkohn ; fi"
   WRITE (UNIT = UnitSh, FMT = "(a)") "echo eigkohn."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "cat eigkohn"
   WRITE (UNIT = UnitSh, FMT = "(a)") "if [ ""$debug"" == """" ] ; then "
   WRITE (UNIT = UnitSh, FMT = "(a)") "   echo modulus from outkohn."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "   egrep ""^ channel momenta:|modulus"" outkohn"
   WRITE (UNIT = UnitSh, FMT = "(a)") "fi"
   WRITE (UNIT = UnitSh, FMT = "(a)") "echo Largest 10 modulus values"
   WRITE (UNIT = UnitSh, FMT = "(a)") "grep modulus outkohn | sort -n -k 8 | tail -10"

   WRITE (UNIT = UnitSh, FMT = "(a)") "mv outkohn outkohn."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv eigkohn eigkohn."//TRIM(ADJUSTL(SymScat))

   WRITE (UNIT = UnitSh, FMT = "(a)") "echo start xcknew a ; date"
   WRITE (UNIT = UnitSh, FMT = "(a)") "cat inkohna >inkohn"
   WRITE (UNIT = UnitSh, FMT = "(a)") "xcknew"
   WRITE (UNIT = UnitSh, FMT = "(a)") "if [ ""$debug"" != """" ] ; then cat outkohn ; fi"
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv outkohn outakohn."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv eigkohn eigakohn."//TRIM(ADJUSTL(SymScat))

   WRITE (UNIT = UnitSh, FMT = "(a)") "echo start xcsolve ; date"
   WRITE (UNIT = UnitSh, FMT = "(a)") "xcsolve"
   WRITE (UNIT = UnitSh, FMT = "(a)") "if [ ""$debug"" != """" ] ; then cat outsolve ; fi"
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv outsolve outsolve."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv insolve insolve."//TRIM(ADJUSTL(SymScat))

   WRITE (UNIT = UnitSh, FMT = "(a)") "}"

   WRITE (UNIT = UnitSh, FMT = "(a)") "function kohnphot {"

   WRITE (UNIT = UnitSh, FMT = "(a)") "echo start xcdipc ; date"
   WRITE (UNIT = UnitSh, FMT = "(a)") "xcdipc"
   WRITE (UNIT = UnitSh, FMT = "(a)") "if [ ""$debug"" != """" ] ; then cat outdipc ; fi"
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv outdipc outdipc."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv indipc indipc."//TRIM(ADJUSTL(SymScat))

   WRITE (UNIT = UnitSh, FMT = "(a)") "echo start xcphot ; date"
   WRITE (UNIT = UnitSh, FMT = "(a)") "xcphot"
   WRITE (UNIT = UnitSh, FMT = "(a)") "if [ ""$debug"" != """" ] ; then cat outphot ; fi"
   WRITE (UNIT = UnitSh, FMT = "(a)") "echo pltphot."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "cat pltphot"
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv outphot outphot."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "cp inphot inphot."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv pltphot pltphot."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv punphot punphot."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "mv dipbf dipbf."//TRIM(ADJUSTL(SymScat))
   WRITE (UNIT = UnitSh, FMT = "(a)") "}"
   CLOSE (UNIT = UnitSh)


! ** write out the geometry information for grid generation program on ingrid

   OPEN (UNIT = UnitTmp, FILE = "ingrid", POSITION = 'REWIND')
   WRITE (UNIT = UnitTmp, FMT = "(i10)") NBuff
   WRITE (UNIT = UnitTmp, FMT = "(3i5)") NCenter, 3, 0
   WRITE (UNIT = UnitTmp, FMT = "(i10)") LMaxK

   Origin = 0.0_XR
   IF (CenterExpand /= 0) THEN
      IF (CenterExpand < 0 .OR. CenterExpand > NCenter) STOP "Bad value for CenterExpand"
      Origin = Coords(:,CenterExpand)
   END IF

   WRITE (UNIT = UnitTmp, FMT = "(3f12.6)") Coords - SPREAD(Origin, 2, NCenter)

   DO i = 1, NCenter
      WRITE (UNIT = UnitTmp, FMT = "('11,1.')")
      WRITE (UNIT = UnitTmp, FMT = "('0.,.5,1.,1.5,2.,3.,5.,10.,30.,50.,100.,150.')")
      WRITE (UNIT = UnitTmp, FMT = "('10,10,10,10,10,20,20,40,40,50,50')")
      WRITE (UNIT = UnitTmp, FMT = "('notheta nophi')")
   END DO

   WRITE (UNIT = UnitTmp, FMT = "('no')")
   CLOSE (UNIT = UnitTmp)

! ** write out the input to the photoionization program in inphot
! ** 
! **    OPEN (UNIT = UnitTmp, FILE = "inphot", POSITION = 'REWIND')
! **    WRITE (UNIT = UnitTmp, FMT = "(3I5)") 0, 0   ! iprt, nww - write out all three transition moments
! **    WRITE (UNIT = UnitTmp, FMT = "(10f12.6)") IPot
! **    CLOSE (UNIT = UnitTmp)

! ** write out list of partial waves to use in inpwlist

   ALLOCATE (npw(Symmetry))
   npw = 0
   DO iSym = 1, Symmetry
      DO L = PWLStart(iSym), LMaxK, CenterSym
         DO M = 0, L
            IF (MOD(M, CnOrder) == PWmMOD(iSym)) THEN
               IF (iCS(iSym) == 0 .OR. M > 0) THEN
                  npw(iSym) = npw(iSym)+1
               END IF
            END IF
         END DO
      END DO
   END DO

   ALLOCATE (PWList(2, MAXVAL(npw), Symmetry))
      
   npw = 0
   DO iSym = 1, Symmetry
      DO L = PWLStart(iSym), LMaxK, CenterSym
         DO M = 0, L
            IF (MOD(M, CnOrder) == PWmMOD(iSym)) THEN
               IF (iCS(iSym) == 0 .OR. M > 0) THEN
                  npw(iSym) = npw(iSym)+1
                  PWList(1, npw(iSym), iSym) = L
                  IF (iCS(iSym) == 0) THEN
                     PWList(2, npw(iSym), iSym) = -M
                  ELSE IF (M > 0) THEN
                     PWList(2, npw(iSym), iSym) = M
                  END IF
               END IF
            END IF
         END DO
      END DO
   END DO

   DO iSym = 1, Symmetry
      OPEN (UNIT = UnitTmp, FILE = "inpwlist."//TRIM(SymCont(iSym)), POSITION = 'REWIND')
      WRITE (UNIT = UnitTmp, FMT = "(i5)") npw(iSym)
      DO i = 1, npw(iSym)
         WRITE (UNIT = UnitTmp, FMT = "(2i5)") PWList(:,i, iSym)
      END DO
      CLOSE (UNIT = UnitTmp)
   END DO

! ** write out the atomic orbital list for the continuum orbital symmetry in inaolist
   WRITE (UNIT = UnitOut, FMT = "('NSym =', i5)") NSym
   DO iSym = 1, Symmetry
      WRITE (UNIT = UnitOut, FMT = "('For symmetry type', i5)") iSym
      WRITE (UNIT = UnitOut, FMT = "('nscont =', i5)") nscont(iSym)
      WRITE (UNIT = UnitOut, FMT = "('FFBFN =', 10i5)") FFBFN(nscont(iSym))
      WRITE (UNIT = UnitOut, FMT = "('ffbfdata =', 10i5)") FFBFData(nscont(iSym))%v
      OPEN (UNIT = UnitTmp, FILE = "inaolist."//TRIM(SymCOnt(iSym)), POSITION = 'REWIND')
      WRITE (UNIT = UnitTmp, FMT = "(i5)") FFBFN(nscont(iSym))
      WRITE (UNIT = UnitTmp, FMT = "(10i5)") FFBFData(nscont(iSym))%v
      CLOSE (UNIT = UnitTmp)
   END DO

! ** write out the molecular orbital list for the continuum orbital symmetry in inmolist

   DO iSym = 1, Symmetry
      OPEN (UNIT = UnitTmp, FILE = "inmolist."//TRIM(SymCont(iSym)), POSITION = 'REWIND')

      WRITE (UNIT = UnitTmp, FMT = "(i5)") NOrbInScat(nscont(iSym))
      WRITE (UNIT = UnitTmp, FMT = "(10i5)") OrbInScat(nscont(iSym))%v
      WRITE (UNIT = UnitTmp, FMT = "(i5)") NOrbInSmall(nscont(iSym))
      WRITE (UNIT = UnitTmp, FMT = "(10i5)") OrbInSmall(nscont(iSym))%v
      CLOSE (UNIT = UnitTmp)
   END DO

   OPEN (UNIT = UnitTmp, FILE = "inmolistcomb", POSITION = 'REWIND')

   WRITE (UNIT = UnitTmp, FMT = "(10i5)") (NOrbInScat(nscont(symroot(ichan))), ichan = 1, nchan)
   DO ichan = 1, nchan
      iSym = symroot(ichan)
      WRITE (UNIT = UnitTmp, FMT = "(10i5)") OrbInScat(nscont(iSym))%v
   END DO
   WRITE (UNIT = UnitTmp, FMT = "(10i5)") (NOrbInSmall(nscont(symroot(ichan))), ichan = 1, nchan)
   DO ichan = 1, nchan
      iSym = symroot(ichan)
      WRITE (UNIT = UnitTmp, FMT = "(10i5)") OrbInSmall(nscont(iSym))%v
   END DO
   CLOSE (UNIT = UnitTmp)
      
! write out data for ffbf in inffbf

   OPEN (UNIT = UnitTmp, FILE = "inffbf", POSITION = 'REWIND')
   WRITE (UNIT = UnitTmp, FMT = "(i10)") NBuff
   WRITE (UNIT = UnitTMp, FMT = "(2i5, /, a)") 0, 1, "tmatrix"
   CLOSE (UNIT = UnitTmp)

   DO ichan = 1, nchan
      iSym = symroot(ichan)
      CALL AppendFile("inffbf", "inpwlist."//TRIM(SymCont(iSym)), UnitTmp2, UnitTmp3)
      CALL AppendFile("inffbf", "inaolist."//TRIM(SymCont(iSym)), UnitTmp2, UnitTmp3)
   END DO
      
! write out data for dipole matrix elements in indipc

   OPEN (UNIT = UnitTmp, FILE = "indipc", POSITION = 'REWIND')
   WRITE (UNIT = UnitTmp, FMT = "(i10)") NBuff
   WRITE (UNIT = UnitTMp, FMT = "(a, /, 2i5)") "coulomb", 0, 1
   CLOSE (UNIT = UnitTmp)

   DO ichan = 1, nchan
      iSym = symroot(ichan)
      CALL AppendFile("indipc", "inpwlist."//TRIM(SymCont(iSym)), UnitTmp2, UnitTmp3)
      CALL AppendFile("indipc", "inmolist."//TRIM(SymCont(iSym)), UnitTmp2, UnitTmp3)
   END DO
! write out data for csolve in insolve

   OPEN (UNIT = UnitTmp, FILE = "insolve", POSITION = 'REWIND')
   WRITE (UNIT = UnitTmp, FMT = "(3i5)") 0, 0, 1
   WRITE (UNIT = UnitTmp, FMT = "(a)") iconq     ! label for contractq
   WRITE (UNIT = UnitTMp, FMT = "(a)") "nexch"
   CLOSE (UNIT = UnitTmp)

! write out data for cknew in inkohna

   OPEN (UNIT = UnitTmp, FILE = "inkohna", POSITION = 'REWIND')
   WRITE (UNIT = UnitTmp, FMT = "(4i5)") 0, 1, 1, 0
! **   WRITE (UNIT = UnitTmp, FMT = "(4i5)") 0, 1, 0, 0  ! no optical potential
   WRITE (UNIT = UnitTMp, FMT = "(a)") "nexch"
   WRITE (UNIT = UnitTMp, FMT = "(a)") "mesa"
   WRITE (UNIT = UnitTMp, FMT = "(a)") "symmetry"
   CLOSE (UNIT = UnitTmp)

   CALL AppendFile("inkohna", "inmolistcomb", UnitTmp2, UnitTmp3)

! write out data for cknew in inkohnb

   OPEN (UNIT = UnitTmp, FILE = "inkohnb", POSITION = 'REWIND')
   WRITE (UNIT = UnitTmp, FMT = "(4i5)") 0, 1, 1, 0
! **    WRITE (UNIT = UnitTmp, FMT = "(4i5)") 0, 1, 0, 0
   WRITE (UNIT = UnitTMp, FMT = "(a)") "nexch"
   WRITE (UNIT = UnitTMp, FMT = "(a)") "bmesa"
   WRITE (UNIT = UnitTMp, FMT = "(a)") "symmetry"
   CLOSE (UNIT = UnitTmp)

   CALL AppendFile("inkohnb", "inmolistcomb", UnitTmp2, UnitTmp3)

! write out data for statnuc in instat

   OPEN (UNIT = UnitTmp, FILE = "instat", POSITION = 'REWIND')
   WRITE (UNIT = UnitTmp, FMT = "(f12.6)") 1.0
   WRITE (UNIT = UnitTmp, FMT = "(i10)") NBuff
   CLOSE (UNIT = UnitTmp)

! write out data for basis in inbas

   OPEN (UNIT = UnitTmp, FILE = "inbas", POSITION = 'REWIND')
   WRITE (UNIT = UnitTmp, FMT = "(i10)") NBuff
   CLOSE (UNIT = UnitTmp)

! write out data for ass in inass

   OPEN (UNIT = UnitTmp, FILE = "inass", POSITION = 'REWIND')
   WRITE (UNIT = UnitTmp, FMT = "(i10, i5, i5)") NBuff, LMaxK, LMaxK
   CLOSE (UNIT = UnitTmp)

! write out data for cbetter1 in incoul

   OPEN (UNIT = UnitTmp, FILE = "incoul", POSITION = 'REWIND')
   WRITE (UNIT = UnitTmp, FMT = "(a)") "coulomb"
   WRITE (UNIT = UnitTmp, FMT = "(f12.6)") 1.0
   WRITE (UNIT = UnitTmp, FMT = "(i5, i5, f12.6, f12.6, i5, f12.6)") LMaxK, 40, 0.0001, 260., 4, 1.0
   CLOSE (UNIT = UnitTmp)

! write out data for ExtractDipMat in inExtractDipMat

   OPEN (UNIT = UnitTmp, FILE = "inExtractDipMat."//TRIM(ADJUSTL(SymScat)), POSITION = 'REWIND')
   WRITE (UNIT = UnitTmp, FMT = "(i5)") NSym
   WRITE (UNIT = UnitTmp, FMT = "(40('''',a5,'''',:,' '))") SymLabelePoly
   WRITE (UNIT = UnitTmp, FMT = "(i5)") nchan
   WRITE (UNIT = UnitTmp, FMT = "(10i5)") (nscont(symroot(ichan)), ichan = 1, nchan)
   WRITE (UNIT = UnitTmp, FMT = "(10i5)") (nstarg(symroot(ichan)), ichan = 1, nchan)
   WRITE (UNIT = UnitTmp, FMT = "(10i5)") nsscat
   CLOSE (UNIT = UnitTmp)

CONTAINS
   SUBROUTINE WriteOccDrt
      IMPLICIT NONE

      INTEGER :: j, i, lastSym, lastGroup

      Line = " "
      LineCnt = 0
      lastSym = 0
      lastGroup = 0
      j = 0

      DO i = 1, nsmall
         IF (GroupList(i) /= lastGroup .OR. OrbSymList(i) /= lastSym) THEN
            IF (j /= 0) THEN
               CALL AddToDRTLine(Line, j, lastGroup, lastSym, LineCnt, LineCntX, UnitTmp)
               j = 0
            END IF
            lastGroup = GroupList(i)
            lastSym = OrbSymList(i)
            j = 1
         ELSE
            j = j+1
         END IF
      END DO
      IF (j /= 0) THEN
         CALL AddToDRTLine(Line, j, lastGroup, lastSym, LineCnt, LineCntX, UnitTmp)
         j = 0
      END IF
   END SUBROUTINE WriteOccDrt
END PROGRAM ExtractData
