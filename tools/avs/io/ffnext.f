*deck @(#)ffnext.f	4.1  7/7/93
      function ffnext(string,pos,start,end)
c***begin prologue     ffnext
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           string, search
c***author             saxe, paul (lanl)
c***source             @(#)ffnext.f	4.1   7/7/93
c***purpose            finds the next token in a string.
c***description
c                      ffnext is a character function used as
c                        found=ffnext(string,pos,start,end)
c                          string   the string to search.
c                          pos      the search begins at pos+1.
c                          start    beginning index of the next token.
c                          end      ending index of the next token.
c
c                      a token is anything delimited by a blank or comma
c                      not within quotes; blanks and commas can therefore
c                      be incorporated within a token by embedding the string
c                      within quotes.
c                      upon return, ffnext denotes what was found:
c                        'eos'          end of string without finding token.
c                        'string'       a character string token.
c                        'integer'      an integer token.
c                        'floating
c                           point'      a floating point token.
c                        'replacement'  a token of the form x=y.
c                        'unknown'      ffnext is confused.
c
c***references
c***routines called    (none)
c***end prologue       ffnext
      implicit integer (a-z)
c
      character*(*) ffnext
      character*(*) string
      character*12 numbs,expon*5
c
      data numbs /'0123456789+-'/, expon /'ddee.'/
      save numbs,expon
c
c
      ffnext='unknown'
c
c     ----- skip blanks and commas -----
c
      do 1 start=pos+1,len(string)
         if((string(start:start).ne.' ')
     $       .and.(string(start:start).ne.',')) goto 2
    1 continue
c
c     ----- nothing but blanks and commas left -----
c
      ffnext='eos'
      start=len(string)
      end=-1
      pos=end
      return
c
c     ----- find end of token, which is a blank or comma not within
c            quotes, or the end of the string
c
    2 continue
      end=start
    3 continue
      newst=end
      junkb=index(string(end:),' ')
      junkc=index(string(end:),',')
      if(junkb.eq.0.or.junkc.eq.0) then
         junk=max(junkb,junkc)
      else
         junk=min(junkb,junkc)
      endif
      if (junk.le.0) then
c
c        ----- no more blanks or commas found, so end of string is
c               end of token --
c
         end=len(string)
      else
         end=end+junk-2
c
c        ----- check if delimiter is within quotes -----
c
         if(end.lt.newst) end=newst
         quote=index(string(newst:end),'"')
         if (quote.gt.0) then
c
c           ----- found a quote, so find matching quote, and look
c                 for a blank again
c
            end=index(string(newst+quote:),'"')
            if (end.le.0) then
               write (6,4) (string(i:i),i=1,len(string))
    4          format (//,' ##### ffnext: unmatched quotes:',/,
     #                (1x,80a1))
               stop
            end if
            end=start+quote+end
            go to 3
         end if
c
      end if
c
c     ----- check for type of object -----
c
      ffnext='string'
      do 11 i=start,end
         if (index(numbs,string(i:i)).eq.0) go to 12
   11 continue
      ffnext='integer'
      pos=end
      return
c
   12 continue
c
c      ----- if starts with d or e consider as string -----
c
       if (index(expon(1:4),string(1:1)).ne.0) go to 14
c
      do 13 i=start,end
         if (index(numbs//expon,string(i:i)).eq.0) go to 14
   13 continue
      ffnext='floating point'
      pos=end
      return
c
   14 continue
      if (index(string(start:end),'=').gt.0) then
         ffnext='replacement'
      end if
c
      pos=end
      return
      end
