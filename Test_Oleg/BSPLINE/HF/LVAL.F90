!=======================================================================
  FUNCTION lval(symbol) 
!=======================================================================
!   Convert symbol to l-value
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER :: lval
    CHARACTER, INTENT(IN) :: symbol
    CHARACTER*22 :: set = 'spdfghiklmnspdfghiklmn'

      locate = index(set,symbol)
      if ( locate .le. 11) then
            lval = locate - 1
         else
            lval = locate - 12
      endif

   END FUNCITON lval
