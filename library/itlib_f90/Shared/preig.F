c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{PREIG: Print Converged Eigenvalues}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck preig.f
c***begin prologue     preig
c***date written       010829   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           residual calculation
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       preig
      subroutine preig(eig,num,begin,ipass)
      implicit integer (a-z)
      real*8 eig
      dimension eig(*)
      common/io/inp, iout
      write(iout,1) ipass
      actual=begin
      do 10 i=1,num
         actual=actual+1
         write(iout,2) actual, eig(i)
 10   continue
      return
 1    format(/,30x,'summary for converged roots on pass = ',i3,
     1       /,17x,'root',22x,'energy')
 2    format(15x,i4,20x,f15.8)
      end

