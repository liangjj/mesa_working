*deck @(#)append.f	1.1  11/30/90
      subroutine append(str,word,pos)
      implicit integer(a-z)
c
c     append the string str to the string word beginning at pos+1.
c     this position is updated before return.
      character*(*) str,word
c
c
      lenstr=len(str)
      word(pos+1:)=str
      pos=pos+lenstr
c
c
      return
      end
