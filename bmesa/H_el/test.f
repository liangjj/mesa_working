c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Non-Iterative Integral Equation without Exchange}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck nonit.f 
c***begin prologue     nonit
c***date written       020303   (yymmdd)
c***revision date               (yymmdd)
c***keywords           integral equation
c***                   
c***author             schneider, b. i.(nsf)
c***source             nonit
c***purpose            noniterative integral equation no exchange
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       nonit
      subroutine nonit(psi,driver,r,r_k,i_k,v_loc,k_int,h,psi_int,n)
      implicit integer (a-z)
      real*8 psi, driver, r, r_k, i_k, v_loc, k_int
      real*8 h, v_psi, h_2, psi_int
      dimension psi(n), driver(n), r(n), r_k(n), i_k(n)
      dimension v_loc(n), k_int(n,2)
      common/io/inp, iout      

c This routine solves the integral equation
c \begin{eqnarray}
c \Psi(r) &=& G(r) + R(r)K^{1}(r) - I(r)K^{2}(r) \\ \nonumber
c K^{1}(r) &=& \int_{0}^{r} dr^{\prime} I(r^{\prime}) V(r^{\prime}) 
c            \Psi(r^{\prime}) \\ \nonumber
c K^{2}(r) &=& \int_{0}^{r} dr^{\prime} R(r^{\prime}) V(r^{\prime}) 
c            \Psi(r^{\prime})
c \end{eqnarray}
c where
c \begin{eqnarray}
c   R(r) = \sqrt{ \frac{2.}{k} } sin(kr) \\ \nonumber
c   I(r) = \sqrt{ \frac{2.}{k} } cos(kr)
c \end{eqnarray}
c Introducing a grid, we can write,
c\begin{eqnarray}
c \Psi_{i} &=& G_{i} + R_{i} K^{1}_{i} - I_{i} K^{2}_{i} \\ \nonumber
c  K^{1}_{i} &=& \int_{0}^{r_{i}} dr^{\prime} I(r^{\prime}) V(r^{\prime}) 
c              \Psi(r^{\prime}) \\ \nonumber
c  K^{2}_{i} &=& \int_{0}^{r_{i}} dr^{\prime} R(r^{\prime}) V(r^{\prime}) 
c              \Psi(r^{\prime}) 
c\end{eqnarray}
c The integrals obey the obvious recursion,
c\begin{eqnarray}
c  K^{1}_{i} &=& K^{1}_{i-1} + \int_{r_{i-1}}^{r_{i}} dr^{\prime} 
c                               I(r^{\prime}) V(r^{\prime}) 
c                              \Psi(r^{\prime}) \\ \nonumber
c  K^{2}_{i} &=& K^{2}_{i-1} + \int_{r_{i-1}}^{r_{i}} dr^{\prime} 
c                               R(r^{\prime}) V(r^{\prime}) 
c                            \Psi(r^{\prime}) 
c\end{eqnarray}
c Using the trapezoidal rule to evaluate the integrals yields,
c\begin{eqnarray}
c  K^{1}_{i} &=& K^{1}_{i-1} + \frac{h}{2} \big [ I_{i} V_{i} \Psi_{i}
c                                        +        I_{i-1} V_{i-1} \Psi_{i-1}
c                                          \big ] \\ \nonumber          
c  K^{2}_{i} &=& K^{2}_{i-1} + \frac{h}{2} \big [ R_{i} V_{i} \Psi_{i}
c                                        +        R_{i-1} V_{i-1} \Psi_{i-1}
c                                          \big ]   
c\end{eqnarray}
c Due to the fact that the kernels are symmetric with respect to the variables
c $r$ and $r^{\prime}$, the contribution of the terms involving the unknown
c function $\Psi_{i}$ cancel.
c So, if we define,
c\begin{eqnarray}
c {\hat K^{1}_{i} } &=& K^{1}_{i-1} + \frac{h}{2} I_{i-1} V_{i-1} 
c                                      \Psi_{i-1} \\ \nonumber
c {\hat K^{2}_{i} } &=& K^{2}_{i-1} + \frac{h}{2} R_{i-1} V_{i-1} 
c                                             \Psi_{i-1}
c\end{eqnarray}
c then the numerical solution of the integral equation does not require 
c the use of the value of the unknown function at the current
c integration point to proceed in the forward direction.  It is this
c property that is crucial to the non-iterative method.
c Initialize the solution, and $K^{1,2}_{i}$ at the first two points.

      h_2=h*.5d0
      k_int(1,1)=0.d0
      k_int(1,2)=0.d0
      psi(1)=driver(1)
      do 10 i=2,n

c        Solution at next point.
c        Calculate $V_{i-1} \Psi_{i-1}$

         v_psi = v_loc(i-1)*psi(i-1) 

c        Perform the partial update for the ${ \hat K^{1,2}_{i} }$
c        integrals needed to compute $\Psi_{1}$.

         k_int(i,1) = k_int(i-1,1) + i_k(i-1)*v_psi*h_2      
         k_int(i,2) = k_int(i-1,2) + r_k(i-1)*v_psi*h_2      

c        compute the solution at the next point.

         psi(i) = driver(i) + r_k(i)*k_int(i,1) - i_k(i)*k_int(i,2)
         v_psi = v_loc(i)*psi(i)
         k_int(i,1) = k_int(i,1) + i_k(i)*v_psi*h_2      
         k_int(i,2) = k_int(i,2) + r_k(i)*v_psi*h_2      
 10   continue   

c     At the end of the integration we store some quantities that are needed
c     later.
c \begin{equation}
c \psi_{int} = K^{1}_{n} 
c \end{equation}

      psi_int = k_int(n,1)      
      return
      end






















