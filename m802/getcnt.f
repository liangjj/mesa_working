*deck @(#)getcnt.f	5.1  11/6/94
      function getcnt()
c
      implicit integer (a-z)
      integer getcnt
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
      if (getcnt.eq.0) getcnt=1
      return
c
    3 getcnt=getcnt*10 + (n-1)
      go to 1
      end
