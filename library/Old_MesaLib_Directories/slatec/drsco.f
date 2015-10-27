*deck drsco
      subroutine drsco (rsav, isav)
c***begin prologue  drsco
c***subsidiary
c***purpose  subsidiary to ddebdf
c***library   slatec
c***type      double precision (rsco-s, drsco-d)
c***author  watts, h. a., (snla)
c***description
c
c   drsco transfers data from arrays to a common block within the
c   integrator package ddebdf.
c
c***see also  ddebdf
c***routines called  (none)
c***common blocks    ddebd1
c***revision history  (yymmdd)
c   820301  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  drsco
c-----------------------------------------------------------------------
c this routine restores from rsav and isav the contents of common
c block ddebd1  , which is used internally in the ddebdf
c package.  this presumes that rsav and isav were loaded by means
c of subroutine dsvco or the equivalent.
c-----------------------------------------------------------------------
c
      integer i, ils, isav, lenils, lenrls
      double precision rls, rsav
      dimension rsav(*),isav(*)
      save lenrls, lenils
      common /ddebd1/ rls(218),ils(33)
      data lenrls /218/, lenils /33/
c
c***first executable statement  drsco
      do 10 i = 1, lenrls
         rls(i) = rsav(i)
   10 continue
      do 20 i = 1, lenils
         ils(i) = isav(i)
   20 continue
      return
c     ----------------------- end of subroutine drsco
c     -----------------------
      end
