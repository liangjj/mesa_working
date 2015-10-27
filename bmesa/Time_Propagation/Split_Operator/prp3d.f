c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Propagate Wavefunction using Split Operator}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle

*deck prp3d.f 
c***begin prologue     prp3d
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, dvr, split-operator
c***                   
c***author             schneider, b. i.(nsf)
c***source             prp3d
c***purpose            split operator propagation.
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       prp3d

      subroutine prp3d(psiout,psiin,psi0,v,work,eig1,eig2,eig3,u1,u2,u3,
     1                 vexp,t,dim,n,nphy)
      implicit integer (a-z)
      real*8 v, eig1, eig2, eig3, u1, u2, u3, t
      real*8 thalf, dum
      complex*16 psiin, psiout, psi0, vexp, work
      dimension nphy(dim)
      dimension psiout(n), psiin(n), psi0(n), v(n)
      dimension eig1(nphy(1)), eig2(nphy(2)), eig3(nphy(3))
      dimension u1(nphy(1),nphy(1)), u2(nphy(2),nphy(2))
      dimension u3(nphy(3),nphy(3))
      dimension vexp(n), work(n)
      common/io/inp, iout     
      thalf = .5d0*t
      call cc2opy(psiin,work,n)

c Step1: Exponentiate the potential and apply to initial state.
c\begin{equation}
c    {\Psi}_{out}(x,y,z,t + {\delta}t) = exp( - i V(x,y,z,t) {\delta}t/2 )
c                                {\Psi}_{in}(x,y,z,t) \nonumber
c\end{equation}
c where, $V(x,y,z,t)$, is the potential.

      call cvexp(vexp,v,thalf,n)
      call cvmul(psiout,vexp,dum,work,n,'complex')   

c Step2: Exponentiate the kinetic energy operator, after transforming to
c the representation which diagonalizes this operator.  This is done for each
c dimension separately.
c

      if(dim.eq.1) then
         call ke1exp(psiout,eig1,u1,t,work,n)
      elseif(dim.eq.2) then
         call ke2exp(psiout,eig1,eig2,u1,u2,t,work,nphy(1),nphy(2))
      elseif(dim.eq.3) then
         call ke3exp(psiout,eig1,eig2,eig3,u1,u2,u3,t,
     1               work,nphy(1),nphy(2),nphy(3))
      else
         call lnkerr('error in ke call')
      endif 

c Step3: Multiply the result by the exponentiated potential.
c        
c

      call cvmul(psiout,vexp,dum,work,n,'complex')

c Step 4: Subtract the right hand side because we are calculting the solution
c         to the inhomgeneous equation and then add back the initial
c         state to get the solution to the Schroedinger equation.

      do 10 i=1,n
         psiout(i) = psiout(i) - psiin(i) + psi0(i)
 10   continue   
c
c Step5: Write the wavefunction to the disk
c

      call iosys ('write real solution to bec',n*2,psiout,0,' ')
      return
      end





















