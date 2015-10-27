*deck dsvco
      subroutine dsvco (rsav, isav)
c***begin prologue  dsvco
c***subsidiary
c***purpose  subsidiary to ddebdf
c***library   slatec
c***type      double precision (svco-s, dsvco-d)
c***author  (unknown)
c***description
c
c   dsvco transfers data from a common block to arrays within the
c   integrator package ddebdf.
c
c***see also  ddebdf
c***routines called  (none)
c***common blocks    ddebd1
c***revision history  (yymmdd)
c   820301  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  dsvco
c-----------------------------------------------------------------------
c this routine stores in rsav and isav the contents of common block
c ddebd1  , which is used internally in the ddebdf package.
c
c rsav = double precision array of length 218 or more.
c isav = integer array of length 33 or more.
c-----------------------------------------------------------------------
      integer i, ils, isav, lenils, lenrls
      double precision rls, rsav
      dimension rsav(*),isav(*)
      save lenrls, lenils
      common /ddebd1/ rls(218),ils(33)
      data lenrls /218/, lenils /33/
c
c***first executable statement  dsvco
      do 10 i = 1, lenrls
         rsav(i) = rls(i)
   10 continue
      do 20 i = 1, lenils
         isav(i) = ils(i)
   20 continue
      return
c     ----------------------- end of subroutine dsvco
c     -----------------------
      end
