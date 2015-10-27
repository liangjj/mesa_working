c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Zfil3}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck zfil3.f
***begin prologue     zfil3
***date written       960723   (yymmdd)
***revision date      yymmdd   (yymmdd)
***keywords           time development
***author             schneider, barry (nsf)
***source             
***purpose            fill array
***references         
***routines called    
***end prologue       zfil3
      subroutine zfil3(psi0,z1,z2,z3,n,m)
      implicit integer (a-z)
      complex*16 psi0, z1, z2, z3
      dimension n(3)
      dimension psi0(m), z1(n(1)), z2(n(2)), z3(n(3))
      common/io/inp, iout
      count=0
      do 10 i=1,n(1)
         do 20 j=1,n(2)
            do 30 k=1,n(3)
               count=count+1
               psi0(count) = z1(i)*z2(j)*z3(k)
 30         continue   
 20      continue
 10   continue   
      return
      end       
