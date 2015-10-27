c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{TSTOVL: Compute Overlap Matrices}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck tstovl.f
      subroutine tstovl(vec,n,nvc)
      implicit integer (a-z)
      real*8 vec, ovl, sdot
      dimension vec(n,nvc)
      common/io/inp, iout
      do 30 i=1,nvc
         do 40 j=1,i
            ovl = sdot(n,vec(1,i),1,vec(1,j),1)
            write(iout,1) i, j, ovl
   40    continue
   30 continue
      return
    1 format(1x,'i = ',i3,1x,'j = ',i3,1x,'overlap = ',e15.8)
      end       
