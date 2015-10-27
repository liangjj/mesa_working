*
*     -----------------------------------------------------------------
*                 l v a l
*     -----------------------------------------------------------------
*
*
      INTEGER FUNCTION lval(symbol)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      CHARACTER*1 symbol
      CHARACTER*22 set
      data         set/'spdfghiklmnspdfghiklmn'/
      locate = index(set,symbol)
      if ( locate .le. 11) then
            lval = locate - 1
         else
            lval = locate - 12
      endif
      return
      end
