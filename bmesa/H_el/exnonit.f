c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Non-Iterative Integral Equation with Exchange}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck exnonit.f 
c***begin prologue     exnonit
c***date written       020303   (yymmdd)
c***revision date               (yymmdd)
c***keywords           integral equation
c***                   
c***author             schneider, b. i.(nsf)
c***source             exnonit
c***purpose            noniterative integral equation with exchange
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       exnonit
      subroutine exnonit(psi,r,psi_1s,driver,r_k,i_k,v_loc,
     1                   k_int,u_int,vpsi_1s,h,psi_int,exch,n)
      implicit integer (a-z)
      real*8 psi, r, psi_1s, driver, r_k, i_k, v_loc, k_int, u_int
      real*8 h, v_psi, vpsi_1s, h_2, psi_int, ex_int
      logical exch
      dimension psi(n), r(n), psi_1s(n), driver(n), r_k(n), i_k(n)
      dimension v_loc(n), k_int(n,2), u_int(n,2), vpsi_1s(n)
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
c Note that $V(r)$ may be non-local and a proper definition is,
c\begin{equation}
c V(r) \Psi(r) = V_{local}(r) \Psi(r) + \int_{0}^{r} dr^{\prime} 
c                          U(r \mid r^{\prime}) \Psi(r^{\prime})
c\end{equation}
c where,
c\begin{eqnarray} 
c \int_{0}^{r} dr^{\prime} U(r \mid r^{\prime}) \Psi(r^{\prime})
c           &=&  \psi_{1s}(r) U^{1}(r) - \frac{ \psi_{1s}(r) } {r} U^{2}(r)
c                          \\ \nonumber
c U^{1}(r) &=& \int_{0}^{r} dr^{\prime}  \frac{ \psi_{1s}(r^{\prime})}{r^{\prime}}
c                                     \Psi(r^{\prime}) \\ \nonumber   
c U^{2}(r) &=& \int_{0}^{r} dr^{\prime} \psi_{1s}(r^{\prime}) 
c                                    \Psi(r^{\prime})
c\end{eqnarray}      
c Introducing a grid, we can write,
c\begin{eqnarray}
c \Psi_{i} &=& G_{i} + R_{i} K^{1}_{i} - I_{i} K^{2}_{i} \\ \nonumber
c  K^{1}_{i} &=& \int_{0}^{r_{i}} dr^{\prime} I(r^{\prime}) V(r^{\prime}) 
c              \Psi(r^{\prime}) \\ \nonumber
c  K^{2}_{i} &=& \int_{0}^{r_{i}} dr^{\prime} R(r^{\prime}) V(r^{\prime}) 
c              \Psi(r^{\prime}) \\ \nonumber
c\end{eqnarray}
c The integrals obey the obvious recursion,
c\begin{eqnarray}
c  K^{1}_{i} &=& K^{1}_{i-1} + \int_{r_{i-1}}^{r_{i}} dr^{\prime} 
c                               I(r^{\prime}) V(r^{\prime}) 
c                              \Psi(r^{\prime}) \\ \nonumber
c  K^{2}_{i} &=& K^{2}_{i-1} + \int_{r_{i-1}}^{r_{i}} dr^{\prime} 
c                               R(r^{\prime}) V(r^{\prime}) 
c                            \Psi(r^{\prime}) \\ \nonumber
c  U^{1}_{i} &=& U^{1}_{i-1} + \int_{r_{i-1}}^{r_{i}} dr^{\prime}                   
c                            \frac{ \psi_{1s}(r^{\prime})}{r^{\prime}}
c                                     \Psi(r^{\prime}) \\ \nonumber
c  U^{2}_{i} &=& U^{2}_{i-1} + \int_{r_{i-1}}^{r_{i}} dr^{\prime} 
c                             \psi_{1s}(r^{\prime}) \Psi(r^{\prime}) 
c\end{eqnarray}
c Using the trapezoidal rule to evaluate the integrals yields,
c\begin{eqnarray}
c  K^{1}_{i} &=& K^{1}_{i-1} + \frac{h}{2} \big [ I_{i} V_{i} \Psi_{i}
c                                        +        I_{i-1} V_{i-1} \Psi_{i-1}
c                                          \big ] \\ \nonumber          
c  K^{2}_{i} &=& K^{2}_{i-1} + \frac{h}{2} \big [ R_{i} V_{i} \Psi_{i}
c                                        +        R_{i-1} V_{i-1} \Psi_{i-1}
c                                          \big ]           \\ \nonumber
c  U^{1}_{i} &=&U^{1}_{i-1} + \frac{h}{2} 
c                            \big [ \frac{ \psi^{i}_{1s} \Psi_{i} } { r_{i} }
c                          +  \frac{ \psi^{i-1}_{1s} \Psi_{i-1} } { r_{i-1} } \big ]
c                                         \\ \nonumber
c  U^{2}_{i} &=& U^{2}_{i-1} + \frac{h}{2} \big [ \psi^{i}_{1s} \Psi_{i}
c                                 +       \psi^{i-1}_{1s} \Psi_{i-1} \big ]
c\end{eqnarray}
c Due to the fact that the kernels are symmetric with respect to the variables
c $r$ and $r^{\prime}$, the contribution of the terms involving the unknown
c function $\Psi_{i}$ cancel.
c So, if we define,
c\begin{eqnarray}
c {\hat K^{1}_{i} } &=& K^{1}_{i-1} + \frac{h}{2} I_{i-1} V_{i-1} 
c                                      \Psi_{i-1} \\ \nonumber
c {\hat K^{2}_{i} } &=& K^{2}_{i-1} + \frac{h}{2} R_{i-1} V_{i-1} 
c                                             \Psi_{i-1} \\ \nonumber
c {\hat U^{1}_{i} } &=& U^{1}_{i-1} + \frac{h}{2} \frac{ \psi^{i-1}_{1s} 
c                                       \Psi_{i-1} } { r^{i-1} } \\ \nonumber
c {\hat U^{2}_{i} } &=& U^{2}_{i-1} + \frac{h}{2} \psi^{i-1}_{1s} 
c                                      \Psi_{i-1} 
c\end{eqnarray}
c then
c\begin{eqnarray}
c \Psi_{i} = G_{i} + R{i} {\hat K^{1}_{i} } - I_{i} {\hat K^{2}_{i}}  \\ \nonumber
c \big [ U \Psi \big ]_{i} = \psi^{i}_{1s} {\hat U^{1}_{i} } 
c                            - \frac{ \psi^{i}_{1s} } { r_{i} } {\hat U^{1}_{i} }
c\end{eqnarray}
c does not require the use of the value of the unknown function at the current
c integration point to proceed in the forward direction.  It is this
c property that is crucial to the non-iterative method.
c Initialize the solution, $K^{1,2}_{i}$ and $U^{1,2}_{i}$ at 
c the first two points.

      h_2=h*.5d0
      k_int(1,1)=0.d0
      k_int(1,2)=0.d0
      u_int(1,1)=0.d0
      u_int(1,2)=0.d0
      psi(1)=driver(1)
      do 10 i=2,n

c        Solution at next point.
c        Calculate $V_{i-1} \Psi_{i-1}$

         v_psi = v_loc(i-1)*psi(i-1) + psi_1s(i-1) * ( u_int(i-1,1)
     1                               - u_int(i-1,2) / r(i-1) )

c        Perform the partial update of the ${\hat K^{1,2}_{i} }$
c        integrals needed to compute $\Psi_{i}$.

         k_int(i,1) = k_int(i-1,1) + i_k(i-1)*v_psi*h_2      
         k_int(i,2) = k_int(i-1,2) + r_k(i-1)*v_psi*h_2      

c        compute the solution at the next point.

         psi(i) = driver(i) + r_k(i)*k_int(i,1) - i_k(i)*k_int(i,2)
         v_psi = v_loc(i)*psi(i)

c        Now do the full update of $U^{1,2}_{i}$

         u_int(i,1) = u_int(i-1,1) + h_2 * ( vpsi_1s(i-1)*psi(i-1)
     1                                   +   vpsi_1s(i)*psi(i) )  
         u_int(i,2) = u_int(i-1,2) + h_2 * ( psi_1s(i-1)*psi(i-1)
     1                                   +   psi_1s(i)*psi(i) )  
         v_psi =  v_psi + psi_1s(i) * u_int(i,1) 
     1                  - vpsi_1s(i) * u_int(i,2)

c        Finish off the update of $K^{1,2}_{i}$

         k_int(i,1) = k_int(i,1) + i_k(i)*v_psi*h_2      
         k_int(i,2) = k_int(i,2) + r_k(i)*v_psi*h_2
 10   continue   

c     At the end of the integration we store, $\langle I \mid V \mid \Psi \rangle$

      psi_int = k_int(n,1)      
      return
      end


















