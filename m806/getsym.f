*deck @(#)getsym.f	5.1  11/6/94
      integer function getsym()
c
      implicit integer (a-z)
c
      character*1 blank,digit,junk,nextc,char
c
      common /tapes/  out,input,drttap
      common /const/  blank, digit(10)
c
      junk=nextc(char)
      do 1 n=1,10
         if (char.eq.digit(n)) go to 3
    1 continue
      call lnkerr('bad symmetry type: '//char)
    3 getsym=n-1-1
      return
      end
