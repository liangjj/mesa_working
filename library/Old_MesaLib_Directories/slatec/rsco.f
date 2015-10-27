*deck rsco
      subroutine rsco (rsav, isav)
c***begin prologue  rsco
c***subsidiary
c***purpose  subsidiary to debdf
c***library   slatec
c***type      single precision (rsco-s, drsco-d)
c***author  watts, h. a., (snla)
c***description
c
c   rsco transfers data from arrays to a common block within the
c   integrator package debdf.
c
c***see also  debdf
c***routines called  (none)
c***common blocks    debdf1
c***revision history  (yymmdd)
c   800901  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  rsco
c
c
c-----------------------------------------------------------------------
c this routine restores from rsav and isav the contents of common
c block debdf1  , which is used internally in the debdf
c package.  this presumes that rsav and isav were loaded by means
c of subroutine svco or the equivalent.
c-----------------------------------------------------------------------
      integer isav, i,      ils, lenils, lenrls
      real rsav, rls
      dimension rsav(*), isav(*)
      common /debdf1/ rls(218), ils(33)
      save lenrls, lenils
      data lenrls/218/, lenils/33/
c
c***first executable statement  rsco
      do 10 i = 1,lenrls
 10     rls(i) = rsav(i)
      do 20 i = 1,lenils
 20     ils(i) = isav(i)
      return
c----------------------- end of subroutine rsco -----------------------
      end
