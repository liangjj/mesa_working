*deck @(#)getcod.f	5.1  11/6/94
      function getcod()
c
      implicit integer (a-z)
      integer getcod
c
      character*1 blank,digit,junk,nextc,multrf,valenc
      character*3 codes,code,words*18
c
      common /tapes/  out,errout,input,drttap
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     #,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     #,               valvir,opensh,multi,speshl
      common /drtchr/ codes(9),words(9),multrf,valenc
      common /const/  blank, digit(10)
c
      junk=nextc(code(1:1))
      junk=nextc(code(2:2))
      junk=nextc(code(3:3))
c
      do 2 i=1,ncodes
         if (code.eq.codes(i)) go to 4
    2 continue
      write (errout,3) code
    3 format (//,' problems with an orbital code:',1x,a3)
      call lnkerr('problems with an orbital code: '//code)
    4 getcod=i
      return
      end
