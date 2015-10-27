*deck xerbla
      subroutine xerbla (srname, info)
c***begin prologue  xerbla
c***subsidiary
c***purpose  error handler for the level 2 and level 3 blas routines.
c***library   slatec
c***category  r3
c***type      all (xerbla-a)
c***keywords  error message
c***author  dongarra, j. j., (anl)
c***description
c
c  purpose
c  =======
c
c  it is called by level 2 and 3 blas routines if an input parameter
c  is invalid.
c
c  parameters
c  ==========
c
c  srname - character*6.
c           on entry, srname specifies the name of the routine which
c           called xerbla.
c
c  info   - integer.
c           on entry, info specifies the position of the invalid
c           parameter in the parameter-list of the calling routine.
c
c***references  (none)
c***routines called  xermsg
c***revision history  (yymmdd)
c   860720  date written
c   910610  routine rewritten to serve as an interface between the
c           level 2 and level 3 blas routines and the slatec error
c           handler xermsg.  (bks)
c***end prologue  xerbla
c
c     ..    scalar arguments ..
      integer            info
      character*6        srname
      character*2        xern1
c
c***first executable statement  xerbla
c
      write (xern1, '(i2)') info
      call xermsg ('slatec', srname, 'on entry to '//srname//
     $             ' parameter number '//xern1//' had an illegal value',
     $             info,1)
c
      return
c
c     end of xerbla.
c
      end
