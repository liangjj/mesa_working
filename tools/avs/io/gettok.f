*deck @(#)gettok.f	4.1  7/7/93
      subroutine gettok(token,string,pos)
c
c***begin prologue     gettok
c***date written       850125   (yymmdd)
c***revision date      860112   (yymmdd)
c***keywords           iosys dependent routines
c
c***author             saxe, paul,    (lanl)
c***source             @(#)gettok.f	4.1   7/7/93
c***purpose            to break the next token (blank delimited) from
c                      a character string.
c***description        #
c
c
c***references
c
c***routines called    lnkerr (mdutil)
c
c   common blocks:     (none)
c
c***end prologue       gettok
c
      implicit integer (a-z)
c
      character*(*) token,string
c
c     ----- if at end of string, return 'eos' -----
c
      if (pos.ge.len(string)) then
         token='eos'
         return
      end if
c
c     ----- skip blanks -----
c
      do 1 i=1,len(string)-pos
         pos=pos+1
         if (string(pos:pos).ne.' ') go to 2
    1 continue
      token='eos'
      return
c
    2 continue
      if (string(pos:pos).ne.'"') then
         start=pos
         pos=pos-1
      else
         start=pos+1
         do 5 i=1,len(string)-pos
            pos=pos+1
            if (string(pos:pos).eq.'"') go to 6
    5    continue
         call lnkerr('gettok: mismatched quotes " ')
    6    continue
         token=string(start:pos-1)
         return
      end if
c
c     ----- find end of token -----
c
      do 3 i=1,len(string)-pos
         pos=pos+1
         if (string(pos:pos).eq.' ') go to 4
    3 continue
      pos=len(string)
    4 continue
      end=min(pos,start+len(token)-1)
      token=string(start:end)
c
      return
      end
