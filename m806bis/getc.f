*deck @(#)getc.f	1.2  7/30/91
      character*1 function getc(char)
c
      implicit integer (a-z)
c
      character*1 char,image
      common /buffer/ pos
      common /buffc/  image(80)
c
      if (pos.ge.80) then
          call getlin
      endif
      pos=pos+1
      char=image(pos)
      getc=char
      return
      end
