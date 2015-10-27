*deck @(#)getcnt.f	1.2  7/30/91
      integer function getcnt()
c
      implicit integer (a-z)
c
      character*1 blank,digit,junk,nextc,char
c
      common /const/  blank, digit(10)
c
      getcnt=0
    1 junk=nextc(char)
      do 2 n=1,10
         if (char.eq.digit(n)) go to 3
    2 continue
      call bkspac
      if (getcnt.eq.0) then
          getcnt=1
      endif
      return
c
    3 getcnt=getcnt*10 + (n-1)
      go to 1
      end
