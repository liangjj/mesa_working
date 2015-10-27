c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Driving Terms}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck drv.f 
c***begin prologue     drv
c***date written       020303   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***                   
c***author             schneider, b. i.(nsf)
c***source             hcoll
c***purpose            driver 
c***                   electron + h atom collisions
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       drv
      subroutine drv(driver,rhs,r_k,i_k,h,n)
c
      implicit integer (a-z)
      real*8 driver, rhs, r_k, i_k
      real*8 h, i_for, i_back, h_2
      dimension driver(n), rhs(n), r_k(n), i_k(n)
      common/io/inp, iout    

c This routine calculates the function
c \begin{eqnarray}
c Driver(r) &=& I(r)K^{1}(r) + R(r)K^{2}(r) \\ \nonumber
c K^{1}(r) &=& \int_{0}^{r} dr^{\prime} R(r^{\prime}) Rhs(r^{\prime}) \\ \nonumber
c K^{2}(r) &=& \int_{r}^{\infty} dr^{\prime} I(r^{\prime}) Rhs(r^{\prime}) 
c            \Psi(r^{\prime})
c \end{eqnarray}
c where $Rhs(r)$ is a known function.

      h_2=.5d0*h
      i_for=0.d0
      driver(1)=0.d0
      do 10 i=2,n
         i_for = i_for + ( r_k(i)*rhs(i) + r_k(i-1)*rhs(i-1) )*h_2
         driver(i) = i_k(i)*i_for
 10   continue   
      i_back=0.d0
      do 20 i=n-1,1,-1
         i_back = i_back + ( i_k(i+1)*rhs(i+1) + i_k(i)*rhs(i) )*h_2
         driver(i) = driver(i) + r_k(i)*i_back
 20   continue   
      return
      end


















