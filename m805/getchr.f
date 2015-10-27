*deck @(#)getchr.f	5.1  11/6/94
      function getchr(char)
c
      implicit integer (a-z)
      character*1 getchr
c
      character*1 char,image
      common /buffer/ pos
      common /buffc/  image(80)
c
      if (pos.ge.80) call getlin
      pos=pos+1
      char=image(pos)
      getchr=char
      return
      end
