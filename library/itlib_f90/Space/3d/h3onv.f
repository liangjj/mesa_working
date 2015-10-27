c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{H1ONV: Matrix Vector Multiply in 3D for DVR Davidson Code}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck h3onv.f
c***begin prologue     h3onv
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
c***end prologue       h3onv
      subroutine h3onv(h1,h2,h3,vecout,vecin,n1,n2,n3,nc,nvc)
      implicit integer (a-z)
      real*8 h1, h2, h3, vecout, vecin
      dimension h1(n1,n1), h2(n2,n2), h3(n3,n3)
      dimension vecin(n3,n2,n1,nc,nvc)
      dimension vecout(n3,n2,n1,nc,nvc)
      common/io/inp, iout
      call apbc(vecout,h3,vecin,n3,n3,n2*n1*nc*nvc)
      do 10 i=1,n1
         do 20 ic=1,nc
            do 30 j=1,nvc
               call apbct(vecout(1,1,i,ic,j),vecin(1,1,i,ic,j),
     1                    h2,n3,n2,n2)
 30         continue   
 20      continue
 10   continue
      do 40 ic=1,nc
         do 50 i=1,nvc
            call apbct(vecout(1,1,1,ic,i),vecin(1,1,1,ic,i),
     1                 h1,n3*n2,n1,n1)
 50      continue   
 40   continue
      return
      end       
