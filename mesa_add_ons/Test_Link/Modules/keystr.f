*deck @(#)keystr.f	5.1  11/6/94
      function keystr(string,key,str,ref)
c
c***begin prologue     keystr
c***date written       860814  (yymmdd)
c***revision date      870131  (yymmdd)
c
c    31 january 1987  pws at lanl
c        removing the declaration of an unused variable 'tmpkey'.
c
c    21 august 1986   modified by pws at lanl
c        adding the ref string and the ability to find truncated keywords
c        if they are present in ref.
c
c***keywords           keyword search
c***author             saxe, paul (lanl)
c***source             @(#)keystr.f	5.1   11/6/94
c***purpose            to return the string associated with a keyword.
c***description
c                character=keystr(string,key,str,ref)
c
c      this function searches for a keyword in a character string and
c      returns an associated substring delimited by , blanks or enclosed
c      in parenthesis. i.e.
c          keyword=substring keyword=....
c          keyword=substring,keyword=....
c          keyword=(sub1,sub2)
c          keyword(sub1,sub2)
c          keyword keyword=....
c          keyword,keyword=....
c
c      there are a couple of special cases: if the keyword is not found, a
c      blank string is return. if the keyword is a standalone, then the
c      keyword itself is returned.
c
c      the string is searched for the entire keyword. if that fails and
c      the keyword is found in the string 'ref', then a search is made
c      on truncated versions of the keyword, until it is found, or it is
c      truncated to such a point that an ambiguity arises with another
c      keyword. in no case will a keyword be truncated to less than 3
c      characters.
c
c      on input:
c               string    character*(*)
c                         the character string to be searched.
c
c               key       character*(*)
c                         the keyword to search for.
c
c               ref       character*(*)
c                         a reference list of keywords for truncation matching.
c                         if you do not wish truncated keywords to match, pass
c                         in a blank string.
c      on return:
c               keystr    character*(*)
c               str       character*(*)
c                         the associated substring.
c
c***references
c***routines called    balpar    (char)
c***end prologue       keystr
c
      implicit integer (a-z)
c
      character*(*) keystr
      character*(*) string,key,str,ref
      character*8 delim
c
      parameter (delim=' ,''()   ')
c
      common /io/     inp,iout
c
c     ----- if the string to search or the key are blank, then
c           it won't be found, so return a blank
c
      if (string.eq.' '.or.key.eq.' ') then
         str=' '
         keystr=' '
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
c-----------------------------------------------------------------------
c  start with the entire key, and see if we can find it in the string.
c  then truncate by one character until either we find it, or come
c  down to an ambiguous length
c-----------------------------------------------------------------------
c
c
c     ----- find the key in the string, if we can -----
c
      start=keyloc(string,key(1:end),junk)
c
      if (start.gt.0) then
c
c        ----- 'next' is the beginning of the rest of the string -----
c
         next=start+end
c
c        ----- if we are at the end of the string, then we have a stand-alone
c
         if (next.gt.len(string)) then
            keystr=key
            str=key
            return
         end if
c
c        ----- if the next character is a delimiter, it is a stand-alone again
c              nb don't consider '(' to be a delimiter in this context
c
         if (index(delim,string(next:next)).ne.0.and.
     #       string(next:next).ne.'(') then
            keystr=key
            str=key
            return
         end if
c
c-------------------------------------------------------------------------
c  at this point the possibilities are:
c    1. '='
c    2. '('
c    3. '=('
c    4. anything else, which means this ain't the right keyword.
c-------------------------------------------------------------------------
c
c        ----- if '=(', shift down one and pick up as '(' -----
c
         if (len(string).gt.next) then
            if (string(next:next+1).eq.'=(') then
               next=next+1
            end if
         end if
c
c        ----- if '(', find matching parenthesis -----
c
         if (string(next:next).eq.'(') then
            if (balpar(string(next:),par).lt.0) then
               write (iout,2) (string(i:i),i=start,min(len(string),
     #                         start+60))
    2          format (//,1x,60a1)
               call lnkerr('balancing parenthesis not found')
            end if
c
            keystr=string(next+1:next+par-2)
            str=string(next+1:next+par-2)
            return
         else if (string(next:next).eq.'=') then
c
c           ----- '=' will give us a simple string, except for quotes
c
            if (len(string).eq.next) then
               call lnkerr('nothing follows keyword= !!')
            end if
            if (string(next+1:next+1).eq.'"') then
               last=index(string(next+2:),'"')
               if (last.eq.0) then
                  write (iout,2) (string(i:i),i=start,min(len(string),
     #                            start+60))
                  call lnkerr('unmatched quotes in keyword="...')
               end if
c
               keystr=string(next+1:next+1+last)
               str=string(next+1:next+1+last)
               return
            else
               do 3 i=next+1,len(string)
                  if (index(delim,string(i:i)).ne.0) then
                     keystr=string(next+1:i-1)
                     str=string(next+1:i-1)
                     return
                  end if
    3          continue
               keystr=string(next+1:)
               str=string(next+1:)
               return
            end if
         end if
      end if
c
c     ----- if the key word is enquoted, we won't check for unambiguous
c              truncations
c
      if (key(1:1).eq.'"') then
         keystr=' '
         str=' '
         return
      end if
c
c     ----- if the key is not in the reference list, the same applies -----
c
      if (keyloc(ref,key(1:end),lockey).eq.0) then
         keystr=' '
         str=' '
         return
      end if
c
c     ----- the minimum disambigous length is arbitrarily 3 except for
c           keys beginning with 'no', which are 5.
c
      if (key(1:2).eq.'no'.or.key(1:2).eq.'no') then
         minlen=5
      else
         minlen=3
      end if
c
c     ----- here begins a truncation loop ------
c
    1 continue
         end=end-1
         if (index(ref(1:lockey-1),key(1:end)).ne.0.or.
     #       index(ref(lockey+1:),key(1:end)).ne.0.or.
     #       end.lt.minlen) then
            keystr=' '
            str=' '
            return
         end if
c
c        ----- find the key in the string, if we can -----
c
         start=keyloc(string,key(1:end),junk)
c
         if (start.gt.0) then
c
c           ----- 'next' is the beginning of the rest of the string -----
c
            next=start+end
c
c           ----- if we are at the end of the string, then we have a stand-alone
c
            if (next.gt.len(string)) then
               keystr=key
               str=key
               return
            end if
c
c           ----- if the next character is a delimiter, it is a stand-alone
c                again nb don't consider '(' to be a delimiter in this context
c
            if (index(delim,string(next:next)).ne.0.and.
     #          string(next:next).ne.'(') then
               keystr=key
               str=key
               return
            end if
c
c-------------------------------------------------------------------------
c  at this point the possibilities are:
c    1. '='
c    2. '('
c    3. '=('
c    4. anything else, which means this ain't the right keyword.
c-------------------------------------------------------------------------
c
c           ----- if '=(', shift down one and pick up as '(' -----
c
            if (len(string).gt.next) then
               if (string(next:next+1).eq.'=(') then
                  next=next+1
               end if
            end if
c
c           ----- if '(', find matching parenthesis -----
c
            if (string(next:next).eq.'(') then
               if (balpar(string(next:),par).lt.0) then
                  write (iout,2) (string(i:i),i=start,min(len(string),
     #                            start+60))
                  call lnkerr('balancing parenthesis not found')
               end if
c
               keystr=string(next:next+par-1)
               str=string(next:next+par-1)
               return
            else if (string(next:next).eq.'=') then
c
c              ----- '=' will give us a simple string, except for quotes
c
               if (len(string).eq.next) then
                  call lnkerr('nothing follows keyword= !!')
               end if
               if (string(next+1:next+1).eq.'"') then
                  last=index(string(next+2:),'"')
                  if (last.eq.0) then
                     write (iout,2) (string(i:i),i=start,
     #                               min(len(string),start+60))
                     call lnkerr('unmatched quotes in keyword="...')
                  end if
c
                  keystr=string(next+1:next+1+last)
                  str=string(next+1:next+1+last)
                  return
               else
                  do 6 i=next+1,len(string)
                     if (index(delim,string(i:i)).ne.0) then
                        keystr=string(next+1:i-1)
                        str=string(next+1:i-1)
                        return
                     end if
    6             continue
                  keystr=string(next+1:)
                  str=string(next+1:)
                  return
               end if
            end if
         end if
      go to 1
c
c
      end
