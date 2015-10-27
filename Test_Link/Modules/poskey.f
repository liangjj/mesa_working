*deck @(#)poskey.f	5.1  11/6/94
      subroutine poskey(string,key,ref,start,end)
c
c***begin prologue     poskey
c***date written       870901  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c
c***keywords           character keyword
c***author             saxe, paul (lanl)
c***source             @(#)poskey.f	5.1   11/6/94
c***purpose            to find the character string associated
c                      with a keyword (or nested keywords)
c***description
c               call poskey(string,key,ref,start,end)
c
c     poskey searches a string for a nested keyword sequence and returns
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
c              ref        character*(*)
c                         the reference string of keywords for truncation
c                         searches. pass in a blank string to disable
c                         searching for truncated keywords.
c
c  on return:
c             start       integer
c                         first location (excluding parenthesis) of key
c                         replacement string.
c
c             end         integer
c                         last location....
c
c***references
c***routines called    nxtkey    (char)
c                      keystr    (char)
c***end prologue       poskey
c
      implicit integer (a-z)
c
      character*(*) string,key,ref
      character*4096 temp2,ref1
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
         call lnkerr('no key passed into poskey !!!')
      end if
c
c     ----- find this key -----
c
      call keypos(string,nkey,ref,start1,end1)
      if (start1.le.0) then
c
c        ----- return default if key not found -----
c
         start=0
         end=-1
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
         call keypos(string(start1:end1),nkey,ref1,start2,end2)
         if (start2.le.0) then
c
c           ----- return default if key not found -----
c
            start=0
            end=-1
            return
         end if
c
c        ---- switch strings and go for next key -----
c
         end1=start1+end2-1
         start1=start1+start2-1
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
      start=start1
      end=end1
c
c
      return
      end
