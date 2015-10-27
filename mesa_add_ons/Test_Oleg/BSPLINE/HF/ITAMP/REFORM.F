*
*     ------------------------------------------------------------------
*               r e f o r m
*     ------------------------------------------------------------------
*
*     convert the free-format str1 to the fixed 5(1x,a3,1x,i4,1x) format
*     for str2
*
      SUBROUTINE reform(str1,str2)
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      INTEGER ouf
      COMMON /INOUT/ iuf,ouf,iscw
      CHARACTER*50 str1,str2,blank
      data                   blank/'   '/
*
    1 i = 0
      str2 = blank
      is = 0
    2 js = index(str1(is+1:),'(')
      if (js .ne. 0) then
         if (js .gt. 5) go to 10
         i = i+5
         str2(i-js+1:i) = str1(is+1:is+js)
         is = is + js
         js = index(str1(is+1:),')')
         if (js .eq. 0 .or. js .gt. 5) go to 10
         i = i+5
         str2(i-js+1:i) = str1(is+1:is+js)
         is = is + js
         go to 2
      end if
      return
   10 write(iscw,*)' error in ',str1,': re-enter'
      read(5,'(a)') str1
      go to 1
      end
