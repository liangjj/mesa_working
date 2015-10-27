*deck @(#)config.f	5.1  11/6/94
      subroutine config(norboc,nbcoc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     norboc is an array of the occupation numbers of the configuration
c     generated, and nbcoc is an array of the basic configuration
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      common /c4/ norb,neli,i,nelu,nexcit,nexctu,ifirst
      dimension norboc(*),nbcoc(*)
      nexctu = 0
      if( ifirst .ne. 0 )   go to 9
      ifirst = 1
      go to 2
c
    9 do 7  k=1,i
      m = nbcoc(k) - norboc(k)
    7 if( m .gt. 0 )   nexctu = nexctu + m
    4 if( norboc(i) .eq. 0 )   go to 8
      norboc(i) = norboc(i) - 1
      nelu = nelu - 1
      if( nbcoc(i) .gt. norboc(i) )   nexctu = nexctu + 1
      if( nexctu .gt. nexcit )   go to 4
c
    2 i = i + 1
      if( i .gt. norb )   go to 1
      norboc(i) = 2
      nelu = nelu + 2
      if( nelu - neli ) 2,5,3
c
    3 norboc(i) = 1
      nelu = nelu - 1
      m = nbcoc(i) - 1
      if( m .gt. 0 )   nexctu = nexctu + m
      if( nexctu .gt. nexcit )   go to 4
c
    5 if( i .eq. norb)   return
      j = i + 1
      mm = 0
      do 6  k=j,norb
      norboc(k) = 0
      m = nbcoc(k) - norboc(k)
      if( m .le. 0 )   go to 6
      nexctu = nexctu + m
      mm = mm + m
    6 continue
      if( nexctu .le. nexcit )   return
      nexctu = nexctu - mm
      go to 4
c
    1 i = i - 1
      nelu = nelu - norboc(i)
      m = nbcoc(i) - norboc(i)
      if( m .gt. 0 )   nexctu = nexctu - m
      i = i - 1
      go to 4
c
    8 m = nbcoc(i) - norboc(i)
      if( m .gt. 0 )   nexctu = nexctu - m
      i = i - 1
      if( i .gt. 0 )   go to 4
      return
      end
