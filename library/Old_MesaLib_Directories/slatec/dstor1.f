*deck dstor1
      subroutine dstor1 (u, yh, v, yp, ntemp, ndisk, ntape)
c***begin prologue  dstor1
c***subsidiary
c***purpose  subsidiary to dbvsup
c***library   slatec
c***type      double precision (stor1-s, dstor1-d)
c***author  watts, h. a., (snla)
c***description
c
c **********************************************************************
c             0 -- storage at output points.
c     ntemp =
c             1 -- temporary storage
c **********************************************************************
c
c***see also  dbvsup
c***routines called  (none)
c***common blocks    dml8sz
c***revision history  (yymmdd)
c   750601  date written
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dstor1
      integer igofx, inhomo, ivp, j, ncomp, nctnf, ndisk, nfc, ntape,
     1     ntemp
      double precision c, u(*), v(*), xsav, yh(*), yp(*)
c
c     ******************************************************************
c
      common /dml8sz/ c,xsav,igofx,inhomo,ivp,ncomp,nfc
c
c      *****************************************************************
c
c     begin block permitting ...exits to 80
c***first executable statement  dstor1
         nctnf = ncomp*nfc
         do 10 j = 1, nctnf
            u(j) = yh(j)
   10    continue
         if (inhomo .eq. 1) go to 30
c
c           zero particular solution
c
c     ......exit
            if (ntemp .eq. 1) go to 80
            do 20 j = 1, ncomp
               v(j) = 0.0d0
   20       continue
         go to 70
   30    continue
c
c           nonzero particular solution
c
            if (ntemp .eq. 0) go to 50
c
               do 40 j = 1, ncomp
                  v(j) = yp(j)
   40          continue
c     .........exit
               go to 80
   50       continue
c
            do 60 j = 1, ncomp
               v(j) = c*yp(j)
   60       continue
   70    continue
c
c        is output information to be written to disk
c
         if (ndisk .eq. 1)
     1      write (ntape) (v(j), j = 1, ncomp),(u(j), j = 1, nctnf)
   80 continue
c
      return
      end
