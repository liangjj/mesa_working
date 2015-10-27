*deck @(#)dot.f	5.1  11/6/94
      function dot(a,b,len)
      implicit integer(a-z)
      real*8 dot,sdot,a(len),b(len)
c
      dot=sdot(len,a,1,b,1)
c
      return
      end
