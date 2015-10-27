*
*   --------------------------------------------------------------------
*               r e o r d
*   --------------------------------------------------------------------
*
*       reorder the list of first appearance so that the FUNCTIONs to be
*   iterated appear last in the list.
*
        SUBROUTINE reord(of, elc, nwf, ierr)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
        CHARACTER*3 of(30), elc
*
        ierr = 1
        call eptr(of, elc, i, *99)
        do 10 j = i, nwf-1
           of(j) = of(j+1)
10      continue
        of(nwf) = elc
        ierr = 0
99      return
        end
