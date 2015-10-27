! ** 67890123 20 890123 30 890123 40 890123 50 890123 60 890123 70 890123 80 890123 90 890123 100 90123 110 90123 120 90
PROGRAM ExtractConfig
! ** rad configuration lists from mesa and count the number of such terms
   USE Extract_mod

   IMPLICIT NONE

   INTEGER (KIND = 4), PARAMETER :: UnitOutConfig = 10, UnitIn = 11, UnitOut = 12, UnitSh = 13
   CHARACTER (LEN = InpLineL) :: Line, symnl2list
   INTEGER (KIND = 4) :: ios
   INTEGER :: Symmetry  ! number of symmetry types
   INTEGER, ALLOCATABLE, DIMENSION(:) :: nrefs, TotalConf
   INTEGER :: iSym, iRef
   INTEGER :: jRef, iCon, iUnq, iTotal
   INTEGER, ALLOCATABLE, DIMENSION(:) :: symncsfs, NOrbInScat
   INTEGER :: ncsfs
   INTEGER :: kref

   OPEN (UNIT = UnitIn, FILE = "inExtractConfig", POSITION = 'REWIND')
   OPEN (UNIT = UnitOut, FILE = "outExtractConfig", POSITION = 'REWIND')
   
   WRITE (UNIT = UnitOut, FMT ="('Program ExtractConfig - get configuration counts from mesa.out')")

   REWIND (UNIT = UnitIn)
   READ (UNIT = UnitIn, FMT = *) Symmetry
   ALLOCATE (nrefs(Symmetry+1))
   READ (UNIT = UnitIn, FMT = *) nrefs
   READ (UNIT = UnitIn, FMT = "(a)") symnl2list
   ALLOCATE (NOrbInScat(Symmetry))
   symnl2list=ADJUSTL(symnl2list)
   READ (UNIT = symnl2list(2:LEN_TRIM(symnl2list)-1), FMT = *) NOrbInScat

   WRITE (UNIT = UnitOut, FMT = "('Number of continuum symmetry types' , i5)") Symmetry
   WRITE (UNIT = UnitOut, FMT = "('Number of numel lines used in each symmetry and for the target/Q space',/, (10i5))")&
        & nrefs
   WRITE (UNIT = UnitOut, FMT = "('Number of basis functions in each continuum', 10i5)") NOrbInScat
   ALLOCATE (TotalConf(Symmetry+1))

   OPEN (UNIT = UnitOutConfig, FILE = "mesa.out", POSITION = 'REWIND')

   DO
      READ (UNIT = UnitOutConfig, IOSTAT = ios, FMT = "(a)") Line
      IF (ios /= 0) THEN
         STOP 'Could not find configuration list'
      END IF
      IF (Line == "         reference   configurations    unique      total    q")  THEN
         EXIT
      END IF
   END DO

   READ (UNIT = UnitOutConfig, FMT = "(a)") Line

   kref = 0
   DO iSym = 1, Symmetry+1
      DO iRef = 1, nrefs(iSym)
         READ (UNIT = UnitOutConfig, FMT = "(a)") Line
         IF (Line(1:32) == "  processing group symmetry info") THEN
            READ (UNIT = UnitOutConfig, FMT = "(a)") Line
            READ (UNIT = UnitOutConfig, FMT = "(a)") Line
         END IF
         READ (UNIT = Line, FMT = *) jRef, iCon, iUnq, iTotal
         kref = kref+1
         IF (jRef /= kRef) THEN
            WRITE (UNIT = UnitOut, FMT = "('Last line read', /, a)") TRIM(Line)
            STOP 'problem reading configuration list'
         END IF
      END DO
      TotalConf(iSym) = iTotal
   END DO
   ncsfs = TotalConf(Symmetry+1)
   WRITE (UNIT = UnitOut, FMT = "('ncsfs', i10)") ncsfs

   DO iSym = Symmetry+1, 2, -1
      TotalConf(iSym) = TotalConf(iSym)-TotalConf(iSym-1)
   END DO
   WRITE (UNIT = UnitOut, FMT = "('Count of each type of reference', /, (10i8))") TotalConf

   ALLOCATE (symncsfs(Symmetry))
   symncsfs = TotalConf(1:Symmetry)/NOrbInScat

   OPEN (UNIT = UnitSh, FILE = "ExtractConfig.sh", POSITION = 'REWIND')
   WRITE (UNIT = UnitSh, FMT = "('#')")

   WRITE (UNIT = UnitSh, FMT = "('symncsfs=""', a , '""  # number for the continuum symmetry basis functions')")&
        & TRIM(MakeList(symncsfs))
   WRITE (UNIT = UnitSh, FMT = "('ncsfs=', a , '  # number of configurations in scattering plus penetration')")&
        & TRIM(ADJUSTL(Int4ToCharL(ncsfs)))

END PROGRAM ExtractConfig
