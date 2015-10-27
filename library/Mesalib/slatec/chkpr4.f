*deck chkpr4
      subroutine chkpr4 (iorder, a, b, m, mbdcnd, c, d, n, nbdcnd, cofx,
     +   idmn, ierror)
c***begin prologue  chkpr4
c***subsidiary
c***purpose  subsidiary to sepx4
c***library   slatec
c***type      single precision (chkpr4-s)
c***author  (unknown)
c***description
c
c     this program checks the input parameters for errors.
c
c***see also  sepx4
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  chkpr4
      external cofx
c***first executable statement  chkpr4
      ierror = 1
      if (a.ge.b .or. c.ge.d) return
c
c     check boundary switches
c
      ierror = 2
      if (mbdcnd.lt.0 .or. mbdcnd.gt.4) return
      ierror = 3
      if (nbdcnd.lt.0 .or. nbdcnd.gt.4) return
c
c     check first dimension in calling routine
c
      ierror = 5
      if (idmn .lt. 7) return
c
c     check m
c
      ierror = 6
      if (m.gt.(idmn-1) .or. m.lt.6) return
c
c     check n
c
      ierror = 7
      if (n .lt. 5) return
c
c     check iorder
c
      ierror = 8
      if (iorder.ne.2 .and. iorder.ne.4) return
c
c     check that equation is elliptic
c
      dlx = (b-a)/m
      do  30 i=2,m
         xi = a+(i-1)*dlx
         call cofx (xi,ai,bi,ci)
      if (ai.gt.0.0) go to 10
      ierror=10
      return
   10 continue
   30 continue
c
c     no error found
c
      ierror = 0
      return
      end
