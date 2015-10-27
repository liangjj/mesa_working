*deck @(#)nxtkey.f	5.1  11/6/94
      function nxtkey(string,pos,key)
c
c***begin prologue     nxtkey
c***date written       860815  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           keyword search
c***author             saxe, paul (lanl)
c***source             @(#)nxtkey.f	5.1   11/6/94
c***purpose            to recover nested keywords.
c***description
c                character=nxtkey(string,pos,key)
c
c     nxtkey will return the next keyword from a nested set of keywords
c     separated by blanks, commas, equals signs, slashes or parenthesis.
c     i.e.   key1=key2=.... or key1 key2 .... or key1,key2,... etc.
c
c   on entry:
c              string    character*(*)
c                        the string containing the keys.
c
c              pos       integer
c                        the location in the string to start searching,
c                        starting with 0.
c
c   on return:
c              nxtkey    character*(*)
c              key       character*(*)
c                        the next key found.
c
c              pos       integer
c                        the ending position of the key found.
c
c   special cases: if there are no more keys, nxtkey and key return a
c                  blank key and pos is set to -1.
c
c***references
c***routines called
c***end prologue      nxtkey
c
      implicit integer (a-z)
c
      character*(*) nxtkey
      character*(*) string,key
      character*8 delim
c
      parameter (delim=' ,=()   ')
c
      common /io/ inp,iout
c
c     ----- check for peculiar parameters -----
c
      if (pos.lt.0) call lnkerr('negative starting position in search')
c
c     ----- if we're at or beyond the ned of the string, we're done -----
c
      if (pos.ge.len(string)) then
         pos=-1
         nxtkey=' '
         key=' '
         return
      end if
c
c     ----- look down the string for a non-delimiter to start key -----
c
      do 1 i=pos+1,len(string)
         if (index(delim,string(i:i)).eq.0) go to 2
    1 continue
c
c     ----- no keywords left if reached here -----
c
      pos=-1
      nxtkey=' '
      key=' '
      return
c
c     ---- found beginning of key, so check for enquoted string -----
c
    2 continue
      begin=i
      if (string(begin:begin).eq.'"') then
         end=index(string(begin+1:),'"')+begin
         if (end.le.begin) then
            write (iout,3) (string(i:i),i=begin,min(len(string),
     #                      begin+60))
    3       format (//,1x,60a1)
            call lnkerr('unmatched quotes in next key')
         end if
c
         pos=end
         nxtkey=string(begin:end)
         key=string(begin:end)
      else
c
c        ----- no quotes, so search for next delimiter -----
c
         do 4 i=begin+1,len(string)
            if (index(delim,string(i:i)).ne.0) go to 5
    4    continue
c
c        ----- key goes to end of string but i is len+1 so next section works
c
    5    continue
         pos=i-1
         key=string(begin:pos)
         nxtkey=string(begin:pos)
      end if
c
c     ----- skip over trailing delimiters -----
c
      do 6 i=pos+2,len(string)
         if (index(delim,string(i:i)).eq.0) then
            pos=i-1
            return
         end if
    6 continue
c
      pos=i-1
c
c
      return
      end
