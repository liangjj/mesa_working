*deck stor1
      subroutine stor1 (u, yh, v, yp, ntemp, ndisk, ntape)
c***begin prologue  stor1
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (stor1-s, dstor1-d)
c***author  watts, h. a., (snla)
c***description
c
c **********************************************************************
c             0 -- storage at output points.
c     ntemp =
c             1 -- temporary storage
c **********************************************************************
c
c***see also  bvsup
c***routines called  (none)
c***common blocks    ml8sz
c***revision history  (yymmdd)
c   750601  date written
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  stor1
      dimension u(*),yh(*),v(*),yp(*)
c
c **********************************************************************
c
      common /ml8sz/ c,xsav,igofx,inhomo,ivp,ncomp,nfc
c
c **********************************************************************
c
c***first executable statement  stor1
      nctnf = ncomp * nfc
      do 10 j = 1,nctnf
   10 u(j) = yh(j)
      if (inhomo .eq. 1)  go to 30
c
c   zero particular solution
c
      if (ntemp .eq. 1)  return
      do 20 j = 1,ncomp
   20 v(j) = 0.
      go to 70
c
c   nonzero particular solution
c
   30 if (ntemp .eq. 0)  go to 50
c
      do 40 j = 1,ncomp
   40 v(j) = yp(j)
      return
c
   50 do 60 j = 1,ncomp
   60 v(j) = c * yp(j)
c
c  is output information to be written to disk
c
   70 if (ndisk .eq. 1)  write (ntape) (v(j),j=1,ncomp),(u(j),j=1,nctnf)
c
      return
      end
