*deck @(#)nextc.f	5.1  11/6/94
      function nextc(char)
c
      implicit character*1 (a-z)
      character*1 nextc
c
      call skipbl
      nextc=getchr(char)
      return
      end
