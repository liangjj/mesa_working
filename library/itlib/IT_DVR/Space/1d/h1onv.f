c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{H1ONV: Matrix Vector Multiply in 1D for DVR Davidson Code}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck h1onv.f
c***begin prologue     h1onv
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            hamiltonian times space vector
c***                   
c***references         
c
c***routines called    
c***end prologue       h1onv
      subroutine h1onv(h1,vecout,vecin,n1,nc,nvc)
      implicit integer (a-z)
      real*8 h1, vecout, vecin
      dimension h1(n1,n1), vecout(n1,nc,nvc), vecin(n1,nc,nvc)
      common/io/inp, iout
      call apbc(vecout,h1,vecin,n1,n1,nc*nvc)
      return
      end       

