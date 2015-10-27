*deck svco
      subroutine svco (rsav, isav)
c***begin prologue  svco
c***subsidiary
c***purpose  subsidiary to debdf
c***library   slatec
c***type      single precision (svco-s, dsvco-d)
c***author  (unknown)
c***description
c
c   svco transfers data from a common block to arrays within the
c   integrator package debdf.
c
c***see also  debdf
c***routines called  (none)
c***common blocks    debdf1
c***revision history  (yymmdd)
c   800901  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  svco
c
c
c-----------------------------------------------------------------------
c this routine stores in rsav and isav the contents of common block
c debdf1  , which is used internally in the debdf package.
c
c rsav = real array of length 218 or more.
c isav = integer array of length 33 or more.
c-----------------------------------------------------------------------
      integer isav, i,      ils, lenils, lenrls
      real rsav, rls
      dimension rsav(*), isav(*)
      common /debdf1/ rls(218), ils(33)
      save lenrls, lenils
      data lenrls/218/, lenils/33/
c
c***first executable statement  svco
      do 10 i = 1,lenrls
 10     rsav(i) = rls(i)
      do 20 i = 1,lenils
 20     isav(i) = ils(i)
      return
c----------------------- end of subroutine svco -----------------------
      end
