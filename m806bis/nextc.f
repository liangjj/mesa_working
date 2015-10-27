*deck @(#)nextc.f	1.2  7/30/91
      character*1 function nextc(char)
c
      implicit character*1 (a-z)
c
      call skipbl
      nextc=getc(char)
      return
      end
