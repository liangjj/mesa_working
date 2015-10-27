*deck @(#)getkey.f	5.1  11/6/94
      function getkey()
c
      implicit integer (a-z)
      character*1 getkey
c
      character*1 blank,digit,junk,nextc,multrf,valenc
      character*3 codes,words*18
c
      common /const/  blank, digit(10)
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     #,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     #,               valvir,opensh,multi,speshl
      common /drtchr/ codes(9),words(9),multrf,valenc
c
      junk=nextc(getkey)
      if (junk.eq.multrf.or.junk.eq.valenc) return
      getkey=blank
      call bkspac
      return
      end
