*deck @(#)getkey.f	1.1  11/30/90
      function getkey(string,icur,keylst,iout)
c***begin prologue     getkey
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           character, search, key
c***author             martin, richard (lanl)
c***source             @(#)getkey.f	1.1   11/30/90
c***purpose            searches for a keyword and it's list.
c***description
c                      getkey is a character function used as:
c                        result=getkey(string,icur,keylst,iout)
c                          string   the string to search.
c                          icur     the search begins at icur+1.
c                          keylst   the list associated with the keyword.
c                          iout     unit number of the print file.
c
c                      a keyword may be specified in three ways.
c                        keyword
c                        keyword=value
c                        keyword="value1,value2,..."
c                      and is separated from other keywords by a blank
c                      or comma.
c                      upon return, the value of getkey contains either the
c                      keyword which was found, or an 'eos' if nothing else
c                      was found in the string.
c
c***references
c***routines called    (ffnext)
c***end prologue       getkey
      implicit integer(a-z)
      character*(*) getkey
      character*(*) string,keylst
      character found*16, ffnext*16, blank*1, quote*1, equals*1
c
      data quote/'"'/, equals/'='/, blank/' '/
c
c
c     look for the next word.
      found=ffnext(string,icur,start,end)
      if(found.eq.'eos') then
c
c        end of line.
         getkey=found
         keylst=blank
         icur=-1
      else if(found.eq.'string') then
c
c        single keyword.
         getkey=string(start:end)
         keylst=blank
         icur=end
      else if(found.eq.'replacement') then
c
c        option(s).
c        find the equals sign.
         poseq=index(string(start:),equals)
         getkey=string(start:start+poseq-2)
c
c        is it a list?
         if(string(start+poseq:start+poseq).eq.quote) then
            posq=index(string(start+poseq+1:),quote)
            keylst=string(start+poseq+1:start+poseq+posq-1)
            icur=end
         else
            keylst=string(start+poseq:end)
            icur=end
         endif
      else
         found=string(start:)
         call lnkerr(' error in keyword syntax.'//found)
      endif
c
c
      return
      end
