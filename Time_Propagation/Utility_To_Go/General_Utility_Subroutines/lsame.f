*deck lsame
      function lsame ( ca, cb )
*     .. scalar arguments ..
      character*1            ca, cb
      logical lsame
*     ..
*
*  purpose
*  =======
*
*  lsame  tests if ca is the same letter as cb regardless of case.
*  cb is assumed to be an upper case letter. lsame returns .true. if
*  ca is either the same as cb or the equivalent lower case letter.
*
*  n.b. this version of the routine is only correct for ascii code.
*       installers must modify the routine for other character-codes.
*
*       for ebcdic systems the constant ioff must be changed to -64.
*       for cdc systems using 6-12 bit representations, the system-
*       specific code in comments must be activated.
*
*  parameters
*  ==========
*
*  ca     - character*1
*  cb     - character*1
*           on entry, ca and cb specify characters to be compared.
*           unchanged on exit.
*
*
*  auxiliary routine for level 2 blas.
*
*  -- written on 20-july-1986
*     richard hanson, sandia national labs.
*     jeremy du croz, nag central office.
*
*     .. parameters ..
      integer                ioff
      parameter            ( ioff=32 )
*     .. intrinsic functions ..
      intrinsic              ichar
*     .. executable statements ..
*
*     test if the characters are equal
*
      lsame = ca .eq. cb
*
*     now test for equivalence
*
      if ( .not.lsame ) then
         lsame = ichar(ca) - ioff .eq. ichar(cb)
      end if
*
      return
*
*  the following comments contain code for cdc systems using 6-12 bit
*  representations.
*
*     .. parameters ..
*     integer                icirfx
*     parameter            ( icirfx=62 )
*     .. scalar arguments ..
*     character*1            cb
*     .. array arguments ..
*     character*1            ca(*)
*     .. local scalars ..
*     integer                ival
*     .. intrinsic functions ..
*     intrinsic              ichar, char
*     .. executable statements ..
*
*     see if the first character in string ca equals string cb.
*
*     lsame = ca(1) .eq. cb .and. ca(1) .ne. char(icirfx)
*
*     if (lsame) return
*
*     the characters are not identical. now check them for equivalence.
*     look for the 'escape' character, circumflex, followed by the
*     letter.
*
*     ival = ichar(ca(2))
*     if (ival.ge.ichar('a') .and. ival.le.ichar('z')) then
*        lsame = ca(1) .eq. char(icirfx) .and. ca(2) .eq. cb
*     end if
*
*     return
*
*     end of lsame.
*
      end


