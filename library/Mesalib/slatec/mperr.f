*deck mperr
      subroutine mperr
c***begin prologue  mperr
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mperr-a)
c***author  (unknown)
c***description
c
c  this routine is called when a fatal error condition is
c  encountered, and after a message has been written on
c  logical unit lun.
c
c***see also  dqdota, dqdoti, mpblas
c***routines called  (none)
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   ??????  modified for use with blas.  blank common changed to named
c           common.  r given dimension 12.
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   930124  increased array size in mpcon for sun -r8.  (rwc)
c***end prologue  mperr
      common /mpcom/ b, t, m, lun, mxr, r(30)
      integer b, t, r
c***first executable statement  mperr
      call xermsg('slatec', 'mperr',
     1   ' *** execution terminated by call to mperr' //
     2   ' in mp version 770217 ***', 1, 2)
c
c at present just stop, but could dump b, t, etc. here.
c action could easily be controlled by a flag in labelled common.
c ansi version uses stop, univac 1108 version uses
c return 0 in order to give a trace-back.
c for debugging purposes it may be useful simply to
c return here.  most mp routines return with result
c zero after calling mperr.
      stop
      end
