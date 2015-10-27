*deck @(#)keyval.f	5.1  11/6/94
      function keyval(string,key,n)
c***begin prologue     keyval
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           string, search, integer, key
c***author             saxe, paul (lanl)
c***source             @(#)keyval.f	5.1   11/6/94
c***purpose            searches for pattern key=integer .
c***description
c                      
c                      searches for a "key=integer" pattern in string.
c                      if the key and value exist, keyval is .true. and
c                      the value is returned in n.  if either the key
c                      or value are missing, keyval is returned .false.
c
c                      keyval is a logical function used as:
c                      test=keyval(string,key,n)
c                        string  the string to search.
c                        substr  the pattern to recognize.
c                        n       value.
c
c***references
c***routines called    ffnext(chr), ctoi(chr)
c***end prologue       keyval
      implicit integer (a-z)
c
      logical keyval
      character*(*) string,key
      character*3 ffnext
c
      keyval=.false.
      loc=index(string,key)
      if (loc.le.0) return
c
      loc=loc+len(key)
      if (loc.gt.len(string)) return
c
      if (string(loc:loc).ne.'=') return
c
      if (ffnext(string,loc,start,end).ne.'int') return
c
      n=ctoi(string(start:end))
      keyval=.true.
c
c
      return
      end
