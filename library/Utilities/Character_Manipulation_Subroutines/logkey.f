*deck @(#)logkey.f	5.1  11/6/94
      function logkey(string,key,default,ref)
c
c***begin prologue     logkey
c***date written       860822  (yymmdd)
c***revision date      870131  (yymmdd)
c
c   31 january 1987  pws at lanl
c
c      fixing a problem with concatenation of key(1:last)//....
c      apparently, since key is assumed length, this is not kosher, so
c      i am transferring 'key' to 'temp', a local 64 character variable.
c      i don't think this is a restriction.
c
c***keywords           logical keyword, standalone keywords
c***author             saxe, paul (lanl)
c***source             @(#)logkey.f	5.1   11/6/94
c***purpose            to search for (negated) stand-alone keywords.
c***description
c                      logkey is a logical function used as:
c               logical=logkey(string,key,default,ref)
c
c     logkey is used to search for standalone keywords, and returns a
c   .true. value if the keyword is found, .false. if the keyword prefixed
c   by 'no' is found, and whatever is passed in as default if the keyword
c   is not found. examples of keywords are 'eigentest' or 'noeigentest'.
c
c   on input:
c              string       character*(*)
c                           the string to search for keywords.
c
c              key          character*(*)
c                           the keyword to search for.
c
c              default      logical
c                           the default value to return if the keyword is
c                           not found.
c
c              ref          character*(*)
c                           the reference string of keywords to use when
c                           searching for an unambiguous truncation of the
c                           keyword.
c
c   on return:
c              logkey       logical
c                           whether the key was found or not, or the default.
c
c***references
c***routines called    keystr
c***end prologue       logkey
c
      implicit integer (a-z)
      logical logkey
c
      character*(*) string,key,ref
      character*64 lastky,temp,nxtkey
      logical default
      logical iskey
c
c     ----- work out what the last subkey is and where it begins -----
c
      pos=0
      if (nxtkey(key,pos,temp).eq.' ') then
c
c        ----- no key passed in so return default -----
c
         logkey=default
         return
      end if
c
      pos1=0
    1 continue
         lastky=temp
         last=pos1
         pos1=pos
      if (nxtkey(key,pos,temp).ne.' ') go to 1
c
c     ----- iskey will return the value of the keyword if it finds
c            the keyword keyword.
c
      if (iskey(string,key,ref)) then
         logkey=.true.
         return
      end if
c
c     ----- negate the keyword and try again -----
c
      if (last.gt.0) then
         temp=key
         if (iskey(string,temp(1:last)//'no'//lastky,ref)) then
            logkey=.false.
            return
         end if
      else
         if (iskey(string,'no'//lastky,ref)) then
            logkey=.false.
            return
         end if
      end if
c
c      ----- haven't found it yet, so default
c
       logkey=default
c
c
       return
       end
