*deck mpchk
      subroutine mpchk (i, j)
c***begin prologue  mpchk
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpchk-a)
c***author  (unknown)
c***description
c
c  checks legality of b, t, m, mxr and lun which should be set
c  in common. the condition on mxr (the dimension of the ep arrays)
c  is that  mxr .ge. (i*t + j)
c
c***see also  dqdota, dqdoti, mpblas
c***routines called  i1mach, mperr
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   ??????  modified for use with blas.  blank common changed to named
c           common.  r given dimension 12.
c   891009  removed unreferenced statement label.  (wrb)
c   891009  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   930124  increased array size in mpcon for sun -r8.  (rwc)
c***end prologue  mpchk
      common /mpcom/ b, t, m, lun, mxr, r(30)
      integer b, t, r
c***first executable statement  mpchk
      lun = i1mach(4)
c now check legality of b, t and m
      if (b.gt.1) go to 40
      write (lun, 30) b
   30 format (' *** b =', i10, ' illegal in call to mpchk,'/
     1 ' perhaps not set before call to an mp routine ***')
      call mperr
   40 if (t.gt.1) go to 60
      write (lun, 50) t
   50 format (' *** t =', i10, ' illegal in call to mpchk,'/
     1 ' perhaps not set before call to an mp routine ***')
      call mperr
   60 if (m.gt.t) go to 80
      write (lun, 70)
   70 format (' *** m .le. t in call to mpchk,'/
     1 ' perhaps not set before call to an mp routine ***')
      call mperr
c 8*b*b-1 should be representable, if not will overflow
c and may become negative, so check for this
   80 ib = 4*b*b - 1
      if ((ib.gt.0).and.((2*ib+1).gt.0)) go to 100
      write (lun, 90)
   90 format (' *** b too large in call to mpchk ***')
      call mperr
c check that space in common is sufficient
  100 mx = i*t + j
      if (mxr.ge.mx) return
c here common is too small, so give error message.
      write (lun, 110) i, j, mx, mxr, t
  110 format (' *** mxr too small or not set to dim(r) before call',
     1 ' to an mp routine *** ' /
     2 ' *** mxr should be at least', i3, '*t +', i4, ' =', i6, '  ***'
     3 / ' *** actually mxr =', i10, ', and t =', i10, '  ***')
      call mperr
      return
      end
