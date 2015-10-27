c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Zfil2}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck zfil2.f
***begin prologue     zfil2
***date written       960723   (yymmdd)
***revision date      yymmdd   (yymmdd)
***keywords           time development
***author             schneider, barry (nsf)
***source             
***purpose            fill array
***references         
***routines called    
***end prologue       zfil2
      subroutine zfil2(psi0,z1,z2,n,m)
      implicit integer (a-z)
      complex*16 psi0, z1, z2
      dimension n(2)
      dimension psi0(m), z1(n(1)), z2(n(2))
      common/io/inp, iout
      count=0
      do 10 i=1,n(1)
         do 20 j=1,n(2)
            count=count+1
            psi0(count) = z1(i)*z2(j)
 20      continue
 10   continue   
      return
      end       


