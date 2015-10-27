*deck @(#)chrkey.f	5.1  11/6/94
      function chrkey(string,key,default,ref)
c
c***begin prologue     chrkey
c***date written       860815  (yymmdd)
c***revision date      860822  (yymmdd)
c
c   22 august 1986  modified by pws at lanl
c        added the 'ref' argument to pass into keystr for searching for
c      truncated keywords.
c
c***keywords           character keyword
c***author             saxe, paul (lanl)
c***source             @(#)chrkey.f	5.1   11/6/94
c***purpose            to find the character string associated
c                      with a keyword (or nested keywords)
c***description
c               character=chrkey(string,key,default,ref)
c
c     chrkey searches a string for a nested keyword sequence and returns
c     the associated character substring. if the keyword is a standalone,
c     then the (last) keyword is returned. if the keyword is not found, the
c     default string is returned.
c
c  on input:
c              string     character*(*)
c                         the string to search.
c
c              key        character*(*)
c                         the key to search for.
c
c              default    character*(*)
c                         the default string to return if the key is not
c                         found.
c
c              ref        character*(*)
c                         the reference string of keywords for truncation
c                         searches. pass in a blank string to disable
c                         searching for truncated keywords.
c
c  on return:
c              chrkey     character*(*)
c                         the associated character string or default string.
c
c***references
c***routines called    nxtkey    (char)
c                      keystr    (char)
c***end prologue       chrkey
c
      implicit integer (a-z)
c
      character*(*) chrkey
      character*(*) string,key,default,ref
      character*4096 temp1,temp2,ref1
      character*60 nkey
      character*1 nxtkey,keystr,junk1
c
c     ----- find the subkeys and there associated strings, one by one
c              the first is different than the rest because use string
c
      pos=0
      if (nxtkey(key,pos,nkey).eq.' ') then
c
c        ----- there is no key to search for -----
c
         call lnkerr('no key passed into chrkey !!!')
      end if
c
c     ----- find this key -----
c
      if (keystr(string,nkey,temp1,ref).eq.' ') then
c
c        ----- return default if key not found -----
c
         chrkey=default
         return
      end if
c
      junk1=keystr(ref,nkey,ref1,' ')
c
c     ----- now loop through keys till no more -----
c
    1 continue
c
c        ----- next key in chain -----
c
         if (nxtkey(key,pos,temp2).eq.' ') go to 2
         nkey=temp2
c
c        ----- and its associated string -----
c
         if (keystr(temp1,nkey,temp2,ref1).eq.' ') then
c
c           ----- return default if key not found -----
c
            chrkey=default
            return
         end if
c
c        ---- switch strings and go for next key -----
c
         temp1=temp2
c
c        ----- limit the reference space -----
c
         junk1=keystr(ref1,nkey,temp2,' ')
         ref1=temp2
      go to 1
c
c     ----- we're at the last keyword and have its associated string -----
c
    2 continue
      if (nkey.ne.temp1) then
         chrkey=temp1
      else
         chrkey=default
      end if
c
c
      return
      end
