c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Exchange Potentials}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck exch.f 
c***begin prologue     exch
c***date written       020303   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***                   
c***author             schneider, b. i.(nsf)
c***source             hcoll
c***purpose            exchange potentials 
c***                   electron + h atom collisions
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       exch
      subroutine exch(ex,psi,psi_1s,r,h,n)
c
      implicit integer (a-z)
      real*8 e, psi, psi_1s, r
      real*8 h, e_for, e_back, h_2
      dimension e(n), psi(n), psi_is(n), r(n)
      common/io/inp, iout    

c This routine calculates the function
c \begin{eqnarray}
c E(r) &=& \psi_{1s}(r) e^{1}(r) + \frac{ \psi_{1s}(r) } {r} e^{2}(r) \\ \nonumber
c e^{1}(r) &=& \int_{0}^{r} dr^{\prime} \frac{ \psi_{1s}(r^{\prime}) 
c                                        \Psi(r^{\prime}) } {r^{\prime}} \\ \nonumber
c e^{2}(r) &=& \int_{r}^{\infty} dr^{\prime} \psi_{1s}(r^{\prime}) \Psi(r^{\prime}) 
c \end{eqnarray}
c where $\Psi(r)$ is a known function.

      h_2=.5d0*h
      e_for=0.d0
      e(1)=0.d0
      do 10 i=2,n
         e_for = e_for + ( psi_1s(i)*psi(i)/r(i) 
     2                +    psi_1s(i-1)*psi(i-1)/r(i-1) )*h_2
         e(i) = psi_1s(i)*e_for
 10   continue   
      e_back=0.d0
      do 20 i=n-1,1,-1
         e_back = e_back + ( psi_1s(i+1)*psi(i+1) + psi_1s(i)*psi(i) )*h_2
         e(i) = e(i) + psi(i)*e_back/r(i)
 20   continue   
      return
      end


















