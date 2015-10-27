*deck @(#)nextc.f	5.1  11/6/94
      character*1 function nextc(char)
c
      implicit character*1 (a-z)
c
      call skipbl
      nextc=getchr(char)
      return
      end
