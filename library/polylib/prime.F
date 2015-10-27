c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Initial Function for Recursion}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle    

*deck prime.f
c***begin prologue     prime
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           scalar product for generalized polynomials
c***author             schneider, barry (nsf)
c***source
c***purpose            
c***                   
c***                   
c***description      
c***                 
c***                                                                       
c***                                                          
c***references         
c
c***routines called    
c***end prologue       prime
      subroutine prime (v,x,n,powr)
      implicit integer (a-z)
      real*8 v, x
      dimension v(n), x(n)
      common/io/inp, iout 
      if(powr.eq.0) then
         call vfill(v,1.d0,n)
      else
         do 10 i=1,n
            v(i)=x(i)**powr
 10      continue   
      endif
      return
      end       
