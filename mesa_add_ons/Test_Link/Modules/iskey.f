*deck @(#)iskey.f	5.1  11/6/94
      function iskey(string,key,ref)
c
c***begin prologue     iskey
c***date written       870721  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c
c***keywords           keywords
c***author             saxe, paul (lanl)
c***source             @(#)iskey.f	5.1   11/6/94
c***purpose            to find if a keyword exists in a string.
c***description
c                      iskey is a logical function used as:
c               logical=iskey(string,key,ref)
c
c     iskey searches a string for a nested keyword sequence and returns
c     whether it finds it.
c
c  on input:
c              string     character*(*)
c                         the string to search.
c
c              key        character*(*)
c                         the key to search for.
c
c              ref        character*(*)
c                         the reference string of keywords for truncation
c                         searches. pass in a blank string to disable
c                         searching for truncated keywords.
c
c  on return:
c              iskey      logical
c                         whether the key exists.
c
c***references
c***routines called    nxtkey    (char)
c                      keystr    (char)
c***end prologue       iskey
c
      implicit integer (a-z)
      logical iskey
c
      character*(*) string,key,ref
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
         call lnkerr('no key passed into iskey !!!')
      end if
c
c     ----- find this key -----
c
      if (keystr(string,nkey,temp1,ref).eq.' ') then
c
c        ----- return default if key not found -----
c
         iskey=.false.
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
            iskey=.false.
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
      iskey=temp1.ne.' '
c
c
      return
      end
