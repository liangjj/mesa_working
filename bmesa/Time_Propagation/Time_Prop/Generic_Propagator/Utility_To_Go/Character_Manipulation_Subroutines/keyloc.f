*deck @(#)keyloc.f	5.1  11/6/94
      function keyloc(string,key,lockey)
c
c***begin prologue     keyloc
c***date written       860821  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           keyword search
c***author             saxe, paul (lanl)
c***source             @(#)keyloc.f	5.1   11/6/94
c***purpose            to return the location in a string of a keyword.
c***description
c                      keyloc is an integer function used as:
c                integer=keyloc(string,key,loc)
c
c      this function searches for a keyword in a character string and
c      returns the location, starting at 1. zero is returned if the string
c      does not contain the key. delimiters are important.
c
c      on input:
c               string    character*(*)
c                         the character string to be searched.
c
c               key       character*(*)
c                         the keyword to search for.
c
c      on return:
c               keyloc    integer
c               lockey    integer
c                         the location of the key.
c
c***references
c***routines called
c***end prologue       keyloc
c
      implicit integer (a-z)
      integer keyloc
c
      character*(*) string,key
      character*8 delim
      character*80 ctemp
c
      external cskipb
c
      parameter (delim=' ,''()=  ')
c
      common /io/ inp,iout
c
c     ----- if the string to search or the key are blank, then
c           it won't be found, so return a blank
c
      if (string.eq.' '.or.key.eq.' ') then
         keyloc=0
         lockey=0
         return
      end if
c
c     ----- search for the end of the key, denoted by a blank. if the
c           key contains blanks it may be enquoted with double quotes (")
c
      if (key(1:1).eq.'"') then
         end=index(key(2:),'"')+1
         if (end.le.1) call lnkerr('unbalanced quotes in key')
      else
         end=index(key,' ')-1
         if (end.lt.0) end=len(key)
      end if
c
c     ----- find the key in the string, if we can -----
c
      begin=1
    1 continue
cnew
         lnb=cskipb(string,' ')
         i=begin
 300     continue
            if (string(i:i+end-1).eq.key(1:end)) then
               start=i-begin+1
               go to 310
            end if
            if (index(delim,string(i:i)).gt.0) then
 301           continue
               if (string(i:i).eq.'(') then
                  offset=balpar(string(i:),junk)
                  if (offset.le.0) then
                     ctemp=string(i:)
                     call lnkerr('no balancing parenthesis found: ...'
     $                    //ctemp)
                  end if
                  i=i+offset
                  go to 300
               else if (string(i:i).eq.'=') then
                  if (string(i:i+1).eq.'=(') then
                     offset=balpar(string(i+1:),junk)
                     if (offset.le.0) then
                        ctemp=string(i:)
                        call lnkerr('no balancing parenthesis found:'//
     $                       ' ...'//ctemp)
                     end if
                     i=i+offset+1
                     go to 300
                  end if
 302              continue
                     i=i+1
                     if (index(delim,string(i:i)).gt.0) go to 301
                  go to 302
               end if
            end if
            i=i+1
         if (i.le.lnb) go to 300
         start=0
 310     continue
cnew
cnew         start=index(string(begin:),key(1:end))
c
c        ----- if it doesn't exist in the string, return blanks -----
c
         if (start.le.0) then
            keyloc=0
            lockey=0
            return
         end if
c
         start=start+begin-1
c
c        ----- 'next' is the beginning of the rest of the string -----
c
         next=start+end
c
c        ----- check for a delimiter in front of the key found to avoid
c               substrings instead of keywords.
c
         if (start.gt.1) then
            if (index(delim,string(start-1:start-1)).eq.0) then
               begin=next
               go to 1
            end if
         end if
c
c        ----- if we are at the end of the string, then we have a stand-alone
c
         if (next.gt.len(string)) then
            keyloc=start
            lockey=start
            return
         end if
c
c        ----- if the next character is a delimiter, it is a stand-alone again
c              nb don't consider '(' to be a delimiter in this context
c
         if (index(delim,string(next:next)).ne.0) then
            keyloc=start
            lockey=start
            return
         end if
c
c        ----- keyword was false one embedded in a larger string -----
c
         begin=next
      go to 1
c
c
      end
