*deck @(#)intkey.f	5.1  11/6/94
      function intkey(string,key,default,ref)
c
c***begin prologue     intkey
c***date written       860815  (yymmdd)
c***revision date      860822  (yymmdd)
c
c   22 august 1986  modified by pws at lanl
c        added 'ref' to the argument list for searching for truncated
c      versions of keywords.
c
c***keywords           integer keyword
c***author             saxe, paul (lanl)
c***source             @(#)intkey.f	5.1   11/6/94
c***purpose            to return an integer associated with a keyword.
c***description
c                      intkey is an integer function used as:
c               integer=intkey(string,key,default,ref)
c
c    intkey searches string for a nested set of keywords in key and returns
c    either the integer value associated with  the keyword if the keyword
c    is found, or the default value if the keyword is not found.
c
c   on input:
c             string      character*(*)
c                         the string to search for the keyword(s).
c
c             key         character*(*)
c                         the keyword(s) to search for.
c
c             default     integer
c                         the default value to return if the keyword is not
c                         found.
c
c             ref         character*(*)
c                         the reference string of keywords for truncation
c                         searches. a blank string disables matching.
c
c   on return:
c             intkey      integer
c                         the integer found, or the default.
c
c***references
c***routines called    chrkey      (char)
c                      ctoi        (char)
c***end prologue       intkey
c
      implicit integer (a-z)
      integer intkey
c
      character*(*) string,key,ref
      character*20 intstr,strtyp,chrkey
c
c     ----- chrkey will return the character representation of the integer,
c           or the keyword for a standalone keyword, or the default we send
c
      intstr=chrkey(string,key,' ',ref)
c
c     ----- blank indicates couldn't find the string -----
c
      if (intstr.eq.' ') then
         intkey=default
         return
      end if
c
c     ----- check for non-integer string, and return default for all such
c           cases
c
      if (strtyp(intstr).ne.'integer') then
         intkey=default
      else
         intkey=ctoi(intstr)
      end if
c
c
      return
      end
