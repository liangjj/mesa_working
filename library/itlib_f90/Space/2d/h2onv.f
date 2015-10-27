c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{H1ONV: Matrix Vector Multiply in 2D for DVR Davidson Code}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck h2onv.f
c***begin prologue     h2onv
c***date written       010828   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            space hamiltonian times space*time vector
c***                   
c***references         
c
c***routines called    
c***end prologue       h2onv
      subroutine h2onv(h1,h2,vecout,vecin,n1,n2,nc,nvc)
      implicit integer (a-z)
      real*8 h1, h2, vecout, vecin
      dimension h1(n1,n1), h2(n2,n2), vecout(n2,n1,nc,nvc)
      dimension vecin(n2,n1,nc,nvc)
      common/io/inp, iout
      call apbc(vecout,h2,vecin,n2,n2,n1*nc*nvc)
      do 10 ic=1,nc
         do 20 i=1,nvc
            call apbct(vecout(1,1,ic,i),vecin(1,1,ic,i),h1,n2,n1,n1)
 20      continue   
 10   continue
      return
      end       
