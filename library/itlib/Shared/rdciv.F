c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{RDCIV: Read in Converged Davidson Vectors}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck rdciv.f
c***begin prologue     rdciv
c***date written       010829   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           read ci vectors
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       rdciv
      subroutine rdciv(vec,code,n,nroot)
      implicit integer (a-z)
      real*8 vec
      character*4 itoc
      character*(*) code
      dimension vec(n,*)
      common/io/inp, iout
      do 10 i=1,nroot
         call iosys('read real "'//code//itoc(i)//'" from rwf',
     1               n,vec(1,i),0,' ')
c         write(iout,1) i, (vec(j,i),j=1,n)
 10   continue
      return
 1    format(/,1x,'ci vector = ',i3,/,5e15.8)
      end

