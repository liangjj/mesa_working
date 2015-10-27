c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Cartesian Kinetic Energy Matrix Elements}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle    

*deck kinxyz.f
c***begin prologue     kinxyz
c***date written       022202   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           kinetic energy
c***author             schneider, barry (nsf)
c***source
c***purpose            generate kinetic energy matrix elements
c***                   in cartesian coordinates.
c***           
c***description
c***references         
c
c***routines called    
c***end prologue       kinxyz
      subroutine kinxyz(kmat,p,ddp,wt,n)
      implicit integer (a-z)
      real*8 kmat, p, ddp, wt
      dimension kmat(n,n), p(n,n), ddp(n,n)
      dimension wt(n)
      common/io/inp, iout 
      do 10 i=1,n
         do 20 j=1,n
            kmat(i,j) = kmat(i,j) + p(i,i)*ddp(i,j)*wt(i)
 20      continue
 10   continue
      return
      end       
