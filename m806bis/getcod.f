*deck @(#)getcod.f	1.2  7/30/91
      integer function getcod()
c
      implicit integer (a-z)
c
      character*1 blank,digit,junk,nextc,multrf,valenc
      character*3 codes,code,words*18
c
      common /tapes/  out,input,drttap
      common /drtcod/ ncodes,dela(9),delb(9),delele(9)
     #,               ntypes,virtul,occupd,valocc,rescor,resvir,frozen
     #,               valvir,opensh,multi,speshl
      common /drtchr/ codes(9),words(9),multrf,valenc
      common /const/  blank, digit(10)
c
      code(1:1)=nextc(junk)
      code(2:2)=nextc(junk)
      code(3:3)=nextc(junk)
c
      if (code(1:3).eq.'typ') then
          code(1:1)=nextc(junk)
          if (nextc(code(2:2)).eq.';') then
              getcod=ctoi(code(1:1))+10
          else if (nextc(code(3:3)).eq.';') then
              getcod=ctoi(code(1:2))+10
          else
              call lnkerr('problems with type code: '//code)
          end if
      else
         do 2 i=1,ncodes
            if (code.eq.codes(i)) go to 4
    2    continue
         call lnkerr(' problems with orbital code: '//code)
    4    getcod=i
      end if
      return
      end
