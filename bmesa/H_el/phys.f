c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Physical Solution}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck phys.f 
c***begin prologue     phys
c***date written       020303   (yymmdd)
c***revision date               (yymmdd)
c***keywords           integral equation
c***                   
c***author             schneider, b. i.(nsf)
c***source             phys
c***purpose            physical solution to scattering integral equation
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       phys
      subroutine phys(psi,psi_00,psi_01,c_00,c_01,n)
      implicit integer (a-z)
      real*8 psi, psi_00, psi_01, c_00, c_01, c_1
      dimension psi(n), psi_00(n), psi_01(n)
      common/io/inp, iout      

c Define the integral operator,
c \begin{equation}
c  G \big( r \mid r^{\prime} \big ) = R(r) I(r^{\prime}) - R(r^{\prime}) I(r)
c\end{equation}
c We wish to find the solution to the equation,
c\begin{eqnarray}
c \Psi(r) &=& D(r) + \int_{0}^{r} dr^{\prime}  
c                     G \big( r \mid r^{\prime} \big ) V(r^{\prime}) 
c                       \Psi(r^{\prime})- R(r) C_{1} \\ \nonumber
c  C_{1} &=& \langle I \mid V \mid \Psi \rangle
c \end{eqnarray}
c The incoming solutions $\Psi^{0}_0$ and $\Psi^{0}_{1}$ satisfy the 
c integral equations,
c \begin{eqnarray}
c \Psi^{0}_{0}(r) &=& R(r) + \int_{0}^{r} dr^{\prime} 
c                     G \big( r \mid r^{\prime} \big ) V(r^{\prime})
c                      \Psi^{0}_{0}(r^{\prime}) \\ \nonumber
c \Psi^{0}_{1}(r) &=& D(r) + \int_{0}^{r} dr^{\prime} 
c                     G \big( r \mid r^{\prime} \big ) V(r^{\prime})
c                     \Psi^{0}_{1}(r^{\prime})
c \end{eqnarray}
c Clearly, since $G$ is a linear operator,
c \begin{equation}
c \Psi(r) = \Psi^{0}_{1}(r) - C_{1} \Psi^{0}_{0}(r)
c \end{equation}
c Defining,
c \begin{eqnarray}
c  C^{0}_{0} &=& \langle I \mid V \mid \Psi^{0}_{0} \rangle \\ \nonumber
c  C^{0}_{1} &=& \langle I \mid V \mid \Psi^{0}_{1} \rangle
c \end{eqnarray}
c gives,
c \begin{equation}
c  C_{1} = \frac{C^{0}_{1}}{ ( 1 + C^{0}_{0} ) }
c \end{equation}

      c_1 = c_01 / ( 1.d0 + c_00 )
      do 10 i=1,n
         psi(i) = psi_01(i) - psi_00(1)*c_1
 10   continue   
      return
      end


















