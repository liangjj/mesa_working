*deck lsame
      logical function lsame (ca, cb)
c***begin prologue  lsame
c***subsidiary
c***purpose  test two characters to determine if they are the same
c            letter, except for case.
c***library   slatec
c***category  r, n3
c***type      logical (lsame-l)
c***keywords  character comparison, level 2 blas, level 3 blas
c***author  hanson, r., (snla)
c           du croz, j., (nag)
c***description
c
c  lsame  tests if ca is the same letter as cb regardless of case.
c  cb is assumed to be an upper case letter. lsame returns .true. if
c  ca is either the same as cb or the equivalent lower case letter.
c
c  n.b. this version of the code is correct for both ascii and ebcdic
c       systems.  installers must modify the routine for other
c       character-codes.
c
c       for cdc systems using 6-12 bit representations, the system-
c       specific code in comments must be activated.
c
c  parameters
c  ==========
c
c  ca     - character*1
c  cb     - character*1
c           on entry, ca and cb specify characters to be compared.
c           unchanged on exit.
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   860720  date written
c   910606  modified to meet slatec prologue standards.  only comment
c           lines were modified.  (bks)
c   910607  modified to handle ascii and ebcdic codes.  (wrb)
c   930201  tests for equality and equivalence combined.  (rwc and wrb)
c***end prologue  lsame
c     .. scalar arguments ..
      character ca*1, cb*1
c     .. local scalars ..
      integer ioff
      logical first
c     .. intrinsic functions ..
      intrinsic ichar
c     .. save statement ..
      save first, ioff
c     .. data statements ..
      data first /.true./
c***first executable statement  lsame
      if (first) ioff = ichar('a') - ichar('a')
c
      first = .false.
c
c     test if the characters are equal or equivalent.
c
      lsame = (ca.eq.cb) .or. (ichar(ca)-ioff.eq.ichar(cb))
c
      return
c
c  the following comments contain code for cdc systems using 6-12 bit
c  representations.
c
c     .. parameters ..
c     integer                icirfx
c     parameter            ( icirfx=62 )
c     .. scalar arguments ..
c     character*1            cb
c     .. array arguments ..
c     character*1            ca(*)
c     .. local scalars ..
c     integer                ival
c     .. intrinsic functions ..
c     intrinsic              ichar, char
c     .. executable statements ..
c     intrinsic              ichar, char
c
c     see if the first character in string ca equals string cb.
c
c     lsame = ca(1) .eq. cb .and. ca(1) .ne. char(icirfx)
c
c     if (lsame) return
c
c     the characters are not identical. now check them for equivalence.
c     look for the 'escape' character, circumflex, followed by the
c     letter.
c
c     ival = ichar(ca(2))
c     if (ival.ge.ichar('a') .and. ival.le.ichar('z')) then
c        lsame = ca(1) .eq. char(icirfx) .and. ca(2) .eq. cb
c     endif
c
c     return
c
c     end of lsame.
c
      end
