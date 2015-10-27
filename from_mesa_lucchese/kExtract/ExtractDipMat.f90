PROGRAM ExtractDipMat
   USE Extract_mod
   IMPLICIT NONE
   INTEGER (KIND = 4), PARAMETER :: Unit_punphot = 11, Unit_dipbf = 12, Unit_inphot = 13, Unit_inExtractDipMat = 14
   INTEGER (KIND = 4), PARAMETER :: Unit_DipMat = 15, UnitOut=16

   INTEGER :: nener, nchan
   INTEGER, ALLOCATABLE, DIMENSION(:) :: nopen
   INTEGER :: idum, i, j, ic, iene, it
   INTEGER, ALLOCATABLE, DIMENSION(:) :: nlm, nscat
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lch, mch, nsch, mePoly
   REAL (KIND = XR) :: eground
   INTEGER :: nmotot
   REAL (KIND = XR), ALLOCATABLE, DIMENSION(:,:) :: hbx, hby, hbz
   REAL (KIND = XR), ALLOCATABLE, DIMENSION(:,:) :: kchan
   INTEGER, ALLOCATABLE, DIMENSION(:) :: iclosed
   COMPLEX (KIND = XR), ALLOCATABLE, DIMENSION(:,:) :: hpvbx, hpvby, hpvbz
   COMPLEX (KIND = XR), ALLOCATABLE, DIMENSION(:,:,:,:) :: tmat

   INTEGER :: ione, itwo, nww, ilm
   REAL (KIND = XR) :: rone
   INTEGER :: iprt
   INTEGER, DIMENSION(3) :: iww
   REAL (KIND = XR), ALLOCATABLE, DIMENSION(:) :: aip
   COMPLEX (KIND = XR), ALLOCATABLE, DIMENSION(:) :: Phase
   INTEGER :: LMaxK
   REAL (KIND = XR) :: ScaleFactor
   CHARACTER (LEN = 2) :: ChanNum
   INTEGER :: NSym
   CHARACTER (LEN = 5), ALLOCATABLE, DIMENSION(:) :: SymLabelePoly
   INTEGER, ALLOCATABLE, DIMENSION(:) :: nscont, nstarg
   INTEGER :: nsscat
   REAL (KIND = XR), DIMENSION(3) :: DipoleOpCoef
   INTEGER :: nchan_inExt
   REAL (KIND = XR) :: PCloseUse
   INTEGER (KIND = 4) :: ios

   OPEN (UNIT = Unit_punphot, FILE="punphot",FORM = "UNFORMATTED", POSITION = "REWIND")
   OPEN (UNIT = Unit_dipbf, FILE="dipbf",FORM = "UNFORMATTED", POSITION = "REWIND")
   OPEN (UNIT = Unit_inphot, FILE="inphot",FORM = "FORMATTED", POSITION = "REWIND")
   OPEN (UNIT = Unit_inExtractDipMat, FILE="inExtractDipMat",FORM = "FORMATTED", POSITION = "REWIND")
   OPEN (UNIT = UnitOut, FILE="outExtractDipMat", FORM="FORMATTED")

   READ (UNIT = Unit_dipbf) nener, nchan, (idum, i = 1, nchan), (idum, i = 1, nchan)
   ALLOCATE (nlm(nchan), nscat(nchan))
   REWIND (UNIT = Unit_dipbf)
   READ (UNIT = Unit_dipbf) nener, nchan, nlm, nscat
   ALLOCATE (lch(MAXVAL(nlm), nchan))
   ALLOCATE (mch(MAXVAL(nlm), nchan),mePoly(MAXVAL(nlm), nchan))
   ALLOCATE (nsch(MAXVAL(nscat), nchan))
   READ (UNIT = Unit_dipbf) ((lch(j, ic), mch(j, ic), j = 1, nlm(ic)), ic = 1, nchan)
   DO ic = 1, nchan
      DO j = 1, nlm(ic)
         mePoly(j, ic) = (mch(j, ic)+1)/2
         IF (2*mePoly(j, ic) == mch(j, ic)) THEN
            mePoly(j, ic) = -mePoly(j, ic)
         END IF
      END DO
   END DO

   READ (UNIT = Unit_dipbf) ((nsch(j, ic), j = 1, nscat(ic)), ic = 1, nchan)
   READ (UNIT = Unit_dipbf) eground
   READ (UNIT = Unit_dipbf) nmotot
   ALLOCATE (hbx(nmotot, nmotot))   
   ALLOCATE (hby(nmotot, nmotot))   
   ALLOCATE (hbz(nmotot, nmotot))
   READ (UNIT = Unit_dipbf) hbx
   READ (UNIT = Unit_dipbf) hby
   READ (UNIT = Unit_dipbf) hbz
   ALLOCATE (kchan(nchan, nener))
   ALLOCATE (iclosed(nchan))

   WRITE (UNIT = UnitOut, FMT = "('Program ExtractDipMat put dipole matrix elements in ePolyScat form')")
   WRITE (UNIT = UnitOut, FMT = "('nchan ', i5)") nchan
   DO ic = 1, nchan
      WRITE (UNIT = UnitOut, FMT = "('Channel ', i5, /, (i5, '  l ', i5, '  m', i5, '  mePoly ', i5))") ic, &
           & (i, lch(i, ic), mch(i, ic), mePoly(i, ic), i = 1, nlm(ic))
   END DO

   READ (UNIT = Unit_inphot, FMT = *) iprt, nww
   IF (nww /= 0) THEN
      READ (UNIT = Unit_inphot, FMT = *) iww(1:nww)
   ELSE
      nww = 3
   END IF
   ALLOCATE (aip(nchan))
   READ (UNIT = Unit_inphot, FMT = *) aip

   READ (UNIT = Unit_inExtractDipMat, FMT = *) NSym
   ALLOCATE (SymLabelePoly(NSym))
   READ (UNIT = Unit_inExtractDipMat, FMT = *) SymLabelePoly
   READ (UNIT = Unit_inExtractDipMat, FMT = *) nchan_inExt
   IF (nchan /= nchan_inExt) THEN 
      STOP "Error - inconsistent nchan values"
   END IF
   ALLOCATE (nscont(nchan), nstarg(nchan))
   READ (UNIT = Unit_inExtractDipMat, FMT = *) nscont
   READ (UNIT = Unit_inExtractDipMat, FMT = *) nstarg
   READ (UNIT = Unit_inExtractDipMat, FMT = *) nsscat
   READ (UNIT = Unit_inExtractDipMat, FMT = *, IOSTAT=ios) PCloseUse
   IF (ios /= 0) THEN
      PCloseUse = PClose
   END IF
   WRITE (UNIT = UnitOut, FMT = "('Ignore blocks of matrix elements if abs of all are less than', e17.8)")&
        & PCloseUse

   ALLOCATE (tmat(3, MAXVAL(nlm), nchan, nener))
   tmat = 0.0_XR
   ALLOCATE (nopen(nener))

   DO iene = 1, nener
      READ (UNIT = Unit_dipbf) kchan(:,iene)
      READ (UNIT = Unit_dipbf) iclosed
      nopen(iene)=0
      DO i = 1, nchan
         IF (iclosed(i) == 0)nopen(iene)=nopen(iene)+1
      END DO
      DO ic = 1, nopen(iene)
         ScaleFactor = SQRT(((aip(ic)/XFEVAU)+0.5_XR*kchan(ic, iene)*kchan(ic, iene))&
              & * XFAUMB*(4.0_XR*PI*PI/(3.0_XR*XFC)))
         ALLOCATE (hpvbx(nlm(ic), nmotot))
         ALLOCATE (hpvby(nlm(ic), nmotot))
         ALLOCATE (hpvbz(nlm(ic), nmotot))
         READ (UNIT = Unit_dipbf) hpvbx
         READ (UNIT = Unit_dipbf) hpvby
         READ (UNIT = Unit_dipbf) hpvbz

         READ (UNIT = Unit_punphot) ione, itwo, nww, ilm, rone, kchan(ic, iene)
         IF (ilm /= nlm(ic)) STOP "Error - Inconsistent ilm"
         IF (nww /= 3) STOP "Error - Expecting all three dipole components"
         READ (UNIT = Unit_punphot) tmat(:,1:nlm(ic), ic, iene)

         LMaxK = MAXVAL(lch(1:nlm(ic), ic))
         ALLOCATE (Phase(0:LMaxK))
         CALL RawPhse(kchan(ic, iene), LMaxK, Phase)
         DO i = 1, nlm(ic)
            tmat(:,i, ic, iene) = -CONJG(tmat(:,i, ic, iene))/Phase(lch(i, ic))
         END DO
         tmat(:, 1:nlm(ic), ic, iene) = tmat(:, 1:nlm(ic), ic, iene)/ScaleFactor

         DO it = 1, 3
            IF (ANY(ABS(tmat(it,1:nlm(ic), ic, iene)) > PCloseUse)) THEN
               WRITE (UNIT = UnitOut, FMT = "('tmat dip =', i5, '  ic =', i5, '  kchan =', e17.8)") &
                    & it, ic, kchan(ic, iene)
               WRITE (UNIT = UnitOut, FMT = "(4e17.8)") tmat(it, 1:nlm(ic), ic, iene)
            END IF
         END DO
         DEALLOCATE (hpvbx, hpvby, hpvbz)
         DEALLOCATE (Phase)
      END DO
   END DO

   DO ic = 1, nchan
      WRITE (UNIT = ChanNum, FMT = "(i2.2)") ic
      OPEN (UNIT = Unit_DipMat, FORM = "FORMATTED", POSITION = "REWIND", FILE = "DipMat"//ChanNum//".dat")
      DO iene = 1, nener
         IF (nopen(iene) >= ic) THEN
            DO it = 1, 3
               IF (ANY(ABS(tmat(it,1:nlm(ic), ic, iene)) > PCloseUse)) THEN
                  WRITE (UNIT = Unit_DipMat, FMT = "(4i10, e16.8, 2i5)") -nscont(ic), 4, 2, nlm(ic), &
                       & 0.5_XR*kchan(ic, iene)*kchan(ic, iene)*XFEVAU, -1, NSym
                  WRITE (UNIT = Unit_DipMat, FMT = "(i10)") NINT(rone)
                  WRITE (UNIT = Unit_DipMat, FMT = "(4i5)") nscont(ic), nstarg(ic), nsscat, -100
                  WRITE (UNIT = Unit_DipMat, FMT = "(3(' ''', a, ''' '))") SymLabelePoly(nscont(ic)),&
                       & SymLabelePoly(nstarg(ic)), SymLabelePoly(nsscat)
                  WRITE (UNIT = Unit_DipMat, FMT = "(4i5)") 1, 1, 1
                  DipoleOpCoef = 0.0_XR
                  DipoleOpCoef(it) = 1.0_XR
                  WRITE (UNIT = Unit_DipMat, FMT = "(4e18.10)") DipoleOpCoef
                  WRITE (UNIT = Unit_DipMat, FMT = "(4e18.10)") 1.0_XR
                  DO i = 1, nlm(ic)
                     WRITE (UNIT = Unit_DipMat, FMT = "(2i5)") lch(i, ic), 1
                     WRITE (UNIT = Unit_DipMat, FMT = "(10i5)") mePoly(i, ic)
                     WRITE (UNIT = Unit_DipMat, FMT = "(4e18.10)") 1.0_XR
                  END DO
                  DO i = 1, nlm(ic)
                     WRITE (UNIT = Unit_DipMat, FMT = "(e20.10, 3e19.10)") tmat(it, i, ic, iene), 0.0_XR, 0.0_XR
                  END DO
               END IF
            END DO
         END IF
      END DO
      CLOSE (UNIT = Unit_DipMat)
   END DO

END PROGRAM ExtractDipMat
