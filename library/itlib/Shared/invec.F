c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{INVEC: Read in CI Vectors for Davidson Code}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck invec.f
c***begin prologue     invec
c***date written       010829   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           read in vectors
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       invec
      subroutine invec(vec,type,n,nvc,prnt)
      implicit integer (a-z)
      real*8 vec
      character*4 itoc
      character*80 title
      character*(*) type
      dimension vec(n,nvc)
      common/io/inp, iout
      do 10 i=1,nvc
         call iosys('read real "'//type//itoc(i)//'" from rwf',
     1               n,vec(1,i),0,' ')
 10   continue
      if(prnt) then
         title='input vectors'
         call prntrm(title,vec,n,nvc,n,nvc,iout)
      endif
      return
      end

