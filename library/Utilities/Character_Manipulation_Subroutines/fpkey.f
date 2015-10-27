*deck @(#)fpkey.f	5.1  11/6/94
      function fpkey(string,key,default,ref)
c
c***begin prologue     fpkey
c***date written       860815  (yymmdd)
c***revision date      860822  (yymmdd)
c
c   22 august 1986  modified by pws at lanl
c        added 'ref' argument for truncated searches to match unambiguous
c      abbreviations of keywords.
c
c***keywords           real*8 keyword, floating point keyword
c***author             saxe, paul (lanl)
c***source             @(#)fpkey.f	5.1   11/6/94
c***purpose            to return a real number associated with a keyword.
c***description
c               real=fpkey(string,key,default,ref)
c
c    fpkey searches string for a nested set of keywords in key and returns
c    either the real value associated with  the keyword if the keyword
c    is found, or the default value if the keyword is not found.
c
c   on input:
c             string      character*(*)
c                         the string to search for the keyword(s).
c
c             key         character*(*)
c                         the keyword(s) to search for.
c
c             default     real
c                         the default value to return if the keyword is not
c                         found.
c
c             ref         character*(*)
c                         the reference string of keywords to use when looking
c                         for an unambiguous abbreviation. pass in a blank
c                         string to disable searching for abbreviated keywords.
c
c   on return:
c             fpkey       real
c                         the integer found, or the default.
c
c***references
c***routines called    chrkey      (char)
c                      ctofp       (char)
c***end prologue       fpkey
c
      implicit integer (a-z)
c
      character*(*) string,key,ref
      character*64 fpstr,strtyp,chrkey
      real*8 default,ctofp,fpkey
c
c     ----- chrkey will return the character representation of the number,
c           or the keyword for a standalone keyword, or the default we send
c
      fpstr=chrkey(string,key,' ',ref)
c
c     ----- blank indicates couldn't find the string -----
c
      if (fpstr.eq.' ') then
         fpkey=default
         return
      end if
c
c     ----- check for non-real string, and return default for all such
c           cases
c
      if (strtyp(fpstr).ne.'real') then
         fpkey=default
      else
         fpkey=ctofp(fpstr)
      end if
c
c
      return
      end
