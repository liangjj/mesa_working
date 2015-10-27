*deck @(#)getsym.f	5.1  11/6/94
      function getsym()
c
      implicit integer (a-z)
      integer getsym
c
      character*1 blank,digit,junk,nextc,char
c
      common /tapes/  out,errout,input,drttap
      common /const/  blank, digit(10)
c
      junk=nextc(char)
      do 1 n=1,10
         if (char.eq.digit(n)) go to 3
    1 continue
      write (errout,2) char
    2 format (//,' funny symmetry type:',1x,a1)
      call lnkerr('bad symmetry type in input: '//char)
    3 getsym=n-1-1
      return
      end
