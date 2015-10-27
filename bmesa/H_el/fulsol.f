c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Construct Physical Solution to Integral Equation}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck fulsol.f 
c***begin prologue     fulsol
c***date written       020303   (yymmdd)
c***revision date               (yymmdd)
c***keywords           integral equation
c***                   
c***author             schneider, b. i.(nsf)
c***source             fulsol
c***purpose            integral equation
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       fulsol
      subroutine fulsol(psi,psi_1s,vpsi_1s,v_lam,e_lam,cmat,c_0,
     1                  ipvt,energy,e_1s,h,nopt,n)
      implicit integer (a-z)
      real*8 psi, psi_1s, vpsi_1s, v_lam, e_lam, cmat, c_0, energy, h
      real*8 e_1s, fac, integ
      dimension psi(n,0:nopt+1), psi_1s(n), vpsi_1s(n)
      dimension v_lam(n,nopt), e_lam(nopt)
      dimension r(n), cmat(nopt+1,nopt+1), c_0(nopt+1), ipvt(*)
      common/io/inp, iout      

c Calculate $\langle \psi_{1s} \mid \frac{1}{r} \mid \Psi_{i} \rangle$,
c $\langle \psi_{1s} \mid  \Psi_{i} \rangle$ and 
c $\langle V_{\lambda} \mid  \Psi_{i} \rangle$ for each value of $i$.
c The solution may be written as,
c \begin{equation}
c \Psi(r) = \Psi_{0}(r) + \big [ \big ( \frac{k^{2}}{2} 
c                                       - \epsilon_{1s} \big )
c           \langle \psi_{1s} \mid \Psi \rangle 
c            - \langle \psi_{1s}\mid \frac{1}{r} \mid \Psi \rangle \big ]
c             \Psi_{1s}(r) - \sum_{\lambda} \langle V_{\lambda} \mid 
c                  \Psi \rangle \frac{\Psi_{\lambda}(r)}{( \epsilon_{1s}
c                                  + \frac{k^2}{2} - \epsilon_{\lambda} )}
c \end{equation}

      fac=energy -e_1s
      ntot=nopt+1
      call rzero(cmat,ntot*ntot)
      call rzero(c_0,ntot)
      c_0(1) = fac * integ(psi_1s,psi(1,0),h,n) 
     1             - 
     2               integ(vpsi_1s,psi(1,0),h,n)
      do 10 i=2,ntot
         c_0(i) = integ(v_lam(1,i-1),psi(1,0),h,n)
 10   continue
      do 20 i=1,ntot   
         cmat(1,i) = fac * integ(psi_1s,psi(1,i),h,n)
     1                   -
     2                     integ(vpsi_1s,psi(1,i),h,n)
         do 30 j=2,ntot
            cmat(j,i) = integ(v_lam(1,j-1),psi(1,i),h,n)
 30      continue   
 20   continue
      call sgefa(cmat,ntot,ntot,ipvt,info)
      call sgesl(cmat,ntot,ntot,ipvt,c_0,0)
      return
      end


















